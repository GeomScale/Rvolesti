#ifndef RANDOM_WALKS_SIMPLEX_BALL_RUN_CHAIN_HPP
#define RANDOM_WALKS_SIMPLEX_BALL_RUN_CHAIN_HPP

#include <stdexcept>

#include <Eigen/Eigen>
#include <cmath>

#include "convex_bodies/simplexintersectball.h"
#include "random_walks/uniform_great_cycle_walk.hpp"
#include "random_walks/uniform_billiard_gcw_walk.hpp"

/// Random walk used by the simplex-ball chain runner.
enum class SimplexBallWalkType
{
    GCW,
    ReGCW
};

namespace simplex_ball_chain_detail
{

template <typename Point, typename MT>
void validate_body_inputs(
    MT const& A,
    Eigen::Matrix<typename Point::FT, Eigen::Dynamic, 1> const& b,
    MT const& V,
    Eigen::Matrix<typename Point::FT, Eigen::Dynamic, 1> const& center)
{
    Eigen::Index const dimension = A.cols();

    if (dimension < 2 ||
        A.rows() <= 0 ||
        b.rows() != A.rows() ||
        V.rows() != dimension ||
        V.cols() <= 0 ||
        center.rows() != dimension)
    {
        throw std::invalid_argument(
            "simplex-ball chain received inconsistent body dimensions");
    }

    if (!A.allFinite() ||
        !b.allFinite() ||
        !V.allFinite() ||
        !center.allFinite())
    {
        throw std::invalid_argument(
            "simplex-ball chain body data must be finite");
    }
}

template <typename Point, typename MT>
SimplexIntersectBall<Point, MT> make_body(
    MT const& A,
    Eigen::Matrix<typename Point::FT, Eigen::Dynamic, 1> const& b,
    MT const& V,
    Eigen::Matrix<typename Point::FT, Eigen::Dynamic, 1> const& center)
{
    typedef typename Point::FT NT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;
    typedef SimplexIntersectBall<Point, MT> Body;

    validate_body_inputs<Point>(A, b, V, center);

    unsigned int const dimension =
        static_cast<unsigned int>(A.cols());

    VT const shifted_b = b - A * center;

    MT shifted_vertices = V;

    for (int vertex = 0; vertex < shifted_vertices.cols(); ++vertex)
    {
        shifted_vertices.col(vertex) -= center;
    }

    return Body(
        dimension,
        A,
        shifted_b,
        shifted_vertices,
        VT::Zero(dimension));
}

template <typename Point, typename Body>
Point make_starting_point(
    Body const& body,
    Eigen::Matrix<typename Point::FT, Eigen::Dynamic, 1> const& start,
    Eigen::Matrix<typename Point::FT, Eigen::Dynamic, 1> const& center)
{
    typedef typename Point::FT NT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;

    if (start.rows() != center.rows() || !start.allFinite())
    {
        throw std::invalid_argument(
            "simplex-ball chain received an invalid starting point");
    }

    VT const shifted_start = start - center;
    NT const squared_radius = shifted_start.squaredNorm();

    if (!std::isfinite(squared_radius) ||
        std::abs(squared_radius - NT(1)) > NT(1e-8))
    {
        throw std::invalid_argument(
            "simplex-ball chain initial point must lie on the unit sphere");
    }

    Point point(shifted_start);

    if (body.is_in(point, NT(1e-10)) == 0)
    {
        throw std::runtime_error(
            "simplex-ball chain initial point is outside the body");
    }

    return point;
}

} // namespace simplex_ball_chain_detail

/// Run either GCW or reflective GCW on one connected component of a
/// simplex-ball intersection, according to walk_type.
///
/// The recentering and sample-collection flow follows
/// R-proj/src/sample_component.cpp; the walk implementations are native
/// volesti C++ policies.
template <
    typename Point,
    typename MT,
    typename RNGType>
MT run_simplex_ball_chain(
    MT const& A,
    Eigen::Matrix<typename Point::FT, Eigen::Dynamic, 1> const& b,
    MT const& V,
    Eigen::Matrix<typename Point::FT, Eigen::Dynamic, 1> const& start,
    unsigned int N,
    unsigned int walk_length,
    Eigen::Matrix<typename Point::FT, Eigen::Dynamic, 1> const& center,
    RNGType& rng,
    unsigned int burnin = 0,
    SimplexBallWalkType walk_type = SimplexBallWalkType::GCW,
    double regcw_tau = 1.0,
    unsigned int max_reflections = 0)
{
    typedef SimplexIntersectBall<Point, MT> Body;

    if (walk_length == 0)
    {
        throw std::invalid_argument(
            "run_simplex_ball_chain: walk length must be positive");
    }

    unsigned int const d = static_cast<unsigned int>(A.cols());

    Body body = simplex_ball_chain_detail::make_body<Point>(
        A,
        b,
        V,
        center);

    Point p = simplex_ball_chain_detail::make_starting_point<Point>(
        body,
        start,
        center);

    // Common burn-in, thinning and sample-collection logic.
    auto collect_samples = [&](auto& walk) -> MT
    {
        // Both walk constructors perform one transition. Count that move as
        // the first burn-in transition, or (when burnin is zero) as the first
        // transition towards the first retained draw. This keeps the public
        // burn-in and thinning counts exact.
        if (burnin > 1)
        {
            walk.apply(body, p, burnin - 1, rng);
        }

        MT samples(d, N);

        if (N == 0)
        {
            return samples;
        }

        unsigned int const first_walk_length =
            burnin == 0
                ? walk_length - 1
                : walk_length;

        if (first_walk_length > 0)
        {
            walk.apply(body, p, first_walk_length, rng);
        }

        samples.col(0) = p.getCoefficients() + center;

        for (unsigned int i = 1; i < N; ++i)
        {
            walk.apply(body, p, walk_length, rng);

            // Move samples back to the original coordinates.
            samples.col(i) = p.getCoefficients() + center;
        }

        return samples;
    };

    if (walk_type == SimplexBallWalkType::GCW)
    {
        typedef GCWalk::template Walk<Body, RNGType> ChainWalk;

        ChainWalk walk(body, p, rng);

        return collect_samples(walk);
    }

    if (walk_type == SimplexBallWalkType::ReGCW)
    {
        if (!std::isfinite(regcw_tau) || regcw_tau <= 0.0)
        {
            throw std::invalid_argument(
                "run_simplex_ball_chain: "
                "ReGCW tau must be positive and finite");
        }

        typedef BilliardGCWalk::template Walk<
            Body,
            RNGType>
            ReflectiveChainWalk;

        if (max_reflections == 0)
        {
            BilliardGCWalk policy(regcw_tau);

            ReflectiveChainWalk walk(
                body,
                p,
                rng,
                policy.param);

            return collect_samples(walk);
        }

        BilliardGCWalk policy(
            regcw_tau,
            max_reflections);

        ReflectiveChainWalk walk(
            body,
            p,
            rng,
            policy.param);

        return collect_samples(walk);
    }

    throw std::invalid_argument(
        "run_simplex_ball_chain: unsupported walk type");
}

#endif
