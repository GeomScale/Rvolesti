// VolEsti (volume computation and sampling library)

// Copyright (c) 2026 Maria Zaza

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef RANDOM_WALKS_SIMPLEX_BALL_CONVERGENCE_HPP
#define RANDOM_WALKS_SIMPLEX_BALL_CONVERGENCE_HPP

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <Eigen/Eigen>

#include "diagnostics/multichain_diagnostics.hpp"
#include "random_walks/simplex_ball_run_chain.hpp"

/// Samples and diagnostics returned by convergence-controlled simplex-ball
/// sampling on one connected component.
template <typename NT, typename VT, typename MT>
struct SimplexBallConvergenceResult
{
    std::vector<MT> chains;
    MT samples;
    MultichainDiagnostics<NT, VT> diagnostics;
    std::vector<MultichainDiagnostics<NT, VT>> diagnostic_history;
    unsigned int diagnostic_checks;
    unsigned int draws_per_chain;
    bool converged;
    std::string warning;
};

namespace simplex_ball_convergence_detail
{

template <typename MT>
MT concatenate_chains(std::vector<MT> const& chains)
{
    Eigen::Index const dimension = chains.front().rows();
    Eigen::Index const draws = chains.front().cols();

    MT samples(
        dimension,
        draws * static_cast<Eigen::Index>(chains.size()));

    for (std::size_t chain = 0; chain < chains.size(); ++chain)
    {
        samples.middleCols(
            static_cast<Eigen::Index>(chain) * draws,
            draws) = chains[chain];
    }

    return samples;
}

template <typename MT>
void append_samples(MT& samples, MT const& new_samples)
{
    if (samples.cols() == 0)
    {
        samples = new_samples;
        return;
    }

    Eigen::Index const previous_draws = samples.cols();

    samples.conservativeResize(
        Eigen::NoChange,
        previous_draws + new_samples.cols());

    samples.rightCols(new_samples.cols()) = new_samples;
}

} // namespace simplex_ball_convergence_detail

/// Run independent chains on one connected simplex-ball component until the
/// Task 4 PSRF/ESS stopping rule is satisfied or max_iterations retained
/// draws per chain have been generated.
///
/// Diagnostics are checked after every check_interval retained draws. Burn-in
/// is applied once per chain and is never included in the diagnostic matrices.
/// All starting points must belong to the same connected component. Supplying
/// points from different components would make Gelman--Rubin diagnostics
/// statistically meaningless.
template <typename Point, typename MT, typename RNGType>
SimplexBallConvergenceResult<
    typename Point::FT,
    Eigen::Matrix<typename Point::FT, Eigen::Dynamic, 1>,
    MT>
run_simplex_ball_chains_until_converged(
    MT const& A,
    Eigen::Matrix<typename Point::FT, Eigen::Dynamic, 1> const& b,
    MT const& V,
    std::vector<
        Eigen::Matrix<typename Point::FT, Eigen::Dynamic, 1>> const&
        starting_points,
    unsigned int check_interval,
    unsigned int max_iterations,
    unsigned int walk_length,
    Eigen::Matrix<typename Point::FT, Eigen::Dynamic, 1> const& center,
    RNGType& rng,
    unsigned int burnin = 0,
    typename Point::FT psrf_target = typename Point::FT(1.1),
    typename Point::FT ess_target = typename Point::FT(100),
    SimplexBallWalkType walk_type = SimplexBallWalkType::GCW,
    double regcw_tau = 1.0,
    unsigned int max_reflections = 0)
{
    typedef typename Point::FT NT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;
    typedef SimplexBallConvergenceResult<NT, VT, MT> Result;

    if (starting_points.size() < 2)
    {
        throw std::invalid_argument(
            "run_simplex_ball_chains_until_converged: "
            "at least two chains are required");
    }

    if (check_interval < 4)
    {
        throw std::invalid_argument(
            "run_simplex_ball_chains_until_converged: "
            "check interval must be at least four");
    }

    if (max_iterations < 4)
    {
        throw std::invalid_argument(
            "run_simplex_ball_chains_until_converged: "
            "maximum iterations must be at least four");
    }

    if (walk_length == 0)
    {
        throw std::invalid_argument(
            "run_simplex_ball_chains_until_converged: "
            "walk length must be positive");
    }

    if (!std::isfinite(psrf_target) || psrf_target <= NT(1))
    {
        throw std::invalid_argument(
            "run_simplex_ball_chains_until_converged: "
            "PSRF target must be finite and above one");
    }

    if (!std::isfinite(ess_target) || ess_target <= NT(0))
    {
        throw std::invalid_argument(
            "run_simplex_ball_chains_until_converged: "
            "ESS target must be positive and finite");
    }

    unsigned int const dimension =
        static_cast<unsigned int>(A.cols());

    if (dimension < 2 ||
        center.rows() != static_cast<Eigen::Index>(dimension))
    {
        throw std::invalid_argument(
            "run_simplex_ball_chains_until_converged: "
            "invalid center dimension");
    }

    for (VT const& starting_point : starting_points)
    {
        if (starting_point.rows() !=
            static_cast<Eigen::Index>(dimension))
        {
            throw std::invalid_argument(
                "run_simplex_ball_chains_until_converged: "
                "invalid starting-point dimension");
        }
    }

    std::size_t const number_of_chains = starting_points.size();

    typedef SimplexIntersectBall<Point, MT> Body;

    Body body = simplex_ball_chain_detail::make_body<Point>(
        A,
        b,
        V,
        center);

    std::vector<Point> points;
    points.reserve(number_of_chains);

    for (VT const& starting_point : starting_points)
    {
        points.push_back(
            simplex_ball_chain_detail::make_starting_point<Point>(
                body,
                starting_point,
                center));
    }

    std::vector<RNGType> chain_rngs;
    chain_rngs.reserve(number_of_chains);

    for (std::size_t chain = 0; chain < number_of_chains; ++chain)
    {
        NT const seed_draw = rng.sample_urdist();

        unsigned int seed =
            static_cast<unsigned int>(
                seed_draw *
                static_cast<NT>(
                    std::numeric_limits<unsigned int>::max()));

        seed ^=
            0x9e3779b9u +
            static_cast<unsigned int>(chain);

        chain_rngs.emplace_back(dimension);
        chain_rngs.back().set_seed(seed);
    }

    std::vector<MT> chains(number_of_chains);
    std::vector<MultichainDiagnostics<NT, VT>> history;
    history.reserve(
        (max_iterations + check_interval - 1) /
        check_interval);

    auto run_sampling_loop = [&](auto& walks) -> Result
    {
        for (std::size_t chain = 0; chain < number_of_chains; ++chain)
        {
            // Construction already performed one transition. Count it as
            // the first burn-in move, rather than silently adding an extra
            // transition before the requested burn-in.
            if (burnin > 1)
            {
                walks[chain].apply(
                    body,
                    points[chain],
                    burnin - 1,
                    chain_rngs[chain]);
            }
        }

        unsigned int draws_per_chain = 0;
        unsigned int diagnostic_checks = 0;
        bool first_draw = true;

        while (draws_per_chain < max_iterations)
        {
            unsigned int const current_batch =
                std::min(
                    check_interval,
                    max_iterations - draws_per_chain);

            for (std::size_t chain = 0;
                 chain < number_of_chains;
                 ++chain)
            {
                MT new_samples(dimension, current_batch);

                for (unsigned int draw = 0;
                     draw < current_batch;
                     ++draw)
                {
                    unsigned int const current_walk_length =
                        first_draw && draw == 0 && burnin == 0
                            ? walk_length - 1
                            : walk_length;

                    if (current_walk_length > 0)
                    {
                        walks[chain].apply(
                            body,
                            points[chain],
                            current_walk_length,
                            chain_rngs[chain]);
                    }

                    new_samples.col(draw) =
                        points[chain].getCoefficients() + center;
                }

                simplex_ball_convergence_detail::append_samples(
                    chains[chain],
                    new_samples);
            }

            first_draw = false;

            draws_per_chain += current_batch;
            ++diagnostic_checks;

            MultichainDiagnostics<NT, VT> diagnostics =
                evaluate_multichain_diagnostics<NT, VT>(
                    chains,
                    psrf_target,
                    ess_target);

            history.push_back(diagnostics);

            if (diagnostics.converged)
            {
                MT samples =
                    simplex_ball_convergence_detail::concatenate_chains(
                        chains);

                return Result{
                    std::move(chains),
                    std::move(samples),
                    diagnostics,
                    std::move(history),
                    diagnostic_checks,
                    draws_per_chain,
                    true,
                    std::string()};
            }
        }

        MultichainDiagnostics<NT, VT> const diagnostics = history.back();

        MT samples =
            simplex_ball_convergence_detail::concatenate_chains(chains);

        return Result{
            std::move(chains),
            std::move(samples),
            diagnostics,
            std::move(history),
            diagnostic_checks,
            draws_per_chain,
            false,
            "maximum iterations reached before PSRF and ESS targets"};
    };

    if (walk_type == SimplexBallWalkType::GCW)
    {
        typedef GCWalk::template Walk<Body, RNGType> ChainWalk;

        std::vector<ChainWalk> walks;
        walks.reserve(number_of_chains);

        for (std::size_t chain = 0;
             chain < number_of_chains;
             ++chain)
        {
            walks.emplace_back(
                body,
                points[chain],
                chain_rngs[chain]);
        }

        return run_sampling_loop(walks);
    }

    if (walk_type == SimplexBallWalkType::ReGCW)
    {
        if (!std::isfinite(regcw_tau) || regcw_tau <= 0.0)
        {
            throw std::invalid_argument(
                "run_simplex_ball_chains_until_converged: "
                "ReGCW tau must be positive and finite");
        }

        typedef BilliardGCWalk::template Walk<Body, RNGType> ChainWalk;

        BilliardGCWalk const policy =
            max_reflections == 0
                ? BilliardGCWalk(regcw_tau)
                : BilliardGCWalk(regcw_tau, max_reflections);

        std::vector<ChainWalk> walks;
        walks.reserve(number_of_chains);

        for (std::size_t chain = 0;
             chain < number_of_chains;
             ++chain)
        {
            walks.emplace_back(
                body,
                points[chain],
                chain_rngs[chain],
                policy.param);
        }

        return run_sampling_loop(walks);
    }

    throw std::invalid_argument(
        "run_simplex_ball_chains_until_converged: unsupported walk type");
}

/// Convenience overload: start all chains from one point, then use independent
/// RNG streams and independent burn-in trajectories.
template <typename Point, typename MT, typename RNGType>
SimplexBallConvergenceResult<
    typename Point::FT,
    Eigen::Matrix<typename Point::FT, Eigen::Dynamic, 1>,
    MT>
run_simplex_ball_chains_until_converged(
    MT const& A,
    Eigen::Matrix<typename Point::FT, Eigen::Dynamic, 1> const& b,
    MT const& V,
    Eigen::Matrix<typename Point::FT, Eigen::Dynamic, 1> const&
        starting_point,
    unsigned int number_of_chains,
    unsigned int check_interval,
    unsigned int max_iterations,
    unsigned int walk_length,
    Eigen::Matrix<typename Point::FT, Eigen::Dynamic, 1> const& center,
    RNGType& rng,
    unsigned int burnin = 0,
    typename Point::FT psrf_target = typename Point::FT(1.1),
    typename Point::FT ess_target = typename Point::FT(100),
    SimplexBallWalkType walk_type = SimplexBallWalkType::GCW,
    double regcw_tau = 1.0,
    unsigned int max_reflections = 0)
{
    typedef typename Point::FT NT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;

    std::vector<VT> starting_points(
        number_of_chains,
        starting_point);

    return run_simplex_ball_chains_until_converged<
        Point,
        MT,
        RNGType>(
            A,
            b,
            V,
            starting_points,
            check_interval,
            max_iterations,
            walk_length,
            center,
            rng,
            burnin,
            psrf_target,
            ess_target,
            walk_type,
            regcw_tau,
            max_reflections);
}

#endif // RANDOM_WALKS_SIMPLEX_BALL_CONVERGENCE_HPP
