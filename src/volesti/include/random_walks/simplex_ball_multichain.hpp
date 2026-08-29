#ifndef RANDOM_WALKS_SIMPLEX_BALL_MULTICHAIN_HPP
#define RANDOM_WALKS_SIMPLEX_BALL_MULTICHAIN_HPP

#include <stdexcept>
#include <vector>
#include <cmath>
#include <cstddef>
#include <limits>
#include <future>

#include <Eigen/Eigen>

#include "random_walks/simplex_ball_run_chain.hpp"

/// Normalize non-negative component volume estimates so that they sum to one.
///
/// This corresponds to:
///     rel_ratios = ratios / sum(ratios)
/// in component_estimation_sampling.R.
template <typename NT>
std::vector<NT> normalize_simplex_ball_component_weights(
    std::vector<NT> const& weights)
{
    if (weights.empty())
    {
        throw std::invalid_argument(
            "normalize_simplex_ball_component_weights: no component weights");
    }

    NT total = NT(0);

    for (NT const& weight : weights)
    {
        if (!std::isfinite(weight) || weight < NT(0))
        {
            throw std::invalid_argument(
                "normalize_simplex_ball_component_weights: "
                "weights must be finite and non-negative");
        }

        total += weight;
    }

    if (!std::isfinite(total) || total <= NT(0))
    {
        throw std::invalid_argument(
            "normalize_simplex_ball_component_weights: "
            "weight sum must be positive and finite");
    }

    std::vector<NT> normalized;
    normalized.reserve(weights.size());

    for (NT const& weight : weights)
    {
        normalized.push_back(weight / total);
    }

    return normalized;
}

/// Draw the component index of every requested sample according to the
/// normalized component volume estimates.
///
/// This corresponds to the runif/cumsum/which logic in
/// component_estimation_sampling.R.
template <typename NT, typename RNGType>
std::vector<std::size_t> sample_simplex_ball_component_assignments(
    std::vector<NT> const& weights,
    unsigned int number_of_samples,
    RNGType& rng)
{
    std::vector<NT> normalized =
        normalize_simplex_ball_component_weights(weights);

    std::vector<NT> cumulative(normalized.size());

    NT running_sum = NT(0);

    for (std::size_t i = 0; i < normalized.size(); ++i)
    {
        running_sum += normalized[i];
        cumulative[i] = running_sum;
    }

    // Avoid a possible floating-point gap below one.
    cumulative.back() = NT(1);

    std::vector<std::size_t> assignments;
    assignments.reserve(number_of_samples);

    for (unsigned int i = 0; i < number_of_samples; ++i)
    {
        NT const u = rng.sample_urdist();
        std::size_t component = 0;

        while (component + 1 < cumulative.size() &&
               u >= cumulative[component])
        {
            ++component;
        }

        assignments.push_back(component);
    }

    return assignments;
}

/// Run one independent Great Cycle Walk chain per selected component and
/// combine the samples according to component volume estimates.
///
/// Unlike the R implementation, which generates a pool of 2*N samples from
/// every component and then resamples from those pools, this implementation
/// first draws the component assignment of each final sample. It then runs
/// every component chain only for the required number of samples.
///
/// raw_component_weights are caller-supplied, finite, non-negative volume
/// estimates with a positive finite sum. They are normalized internally.
/// Components assigned no output samples are skipped and no chain is started
/// for them.
template <
    typename Point,
    typename MT,
    typename RNGType>
MT run_simplex_ball_multichain(
    MT const& A,
    Eigen::Matrix<typename Point::FT, Eigen::Dynamic, 1> const& b,
    MT const& V,
    std::vector<
        Eigen::Matrix<typename Point::FT, Eigen::Dynamic, 1>> const&
        starting_points,
    std::vector<typename Point::FT> const& raw_component_weights,
    unsigned int number_of_samples,
    unsigned int walk_length,
    Eigen::Matrix<typename Point::FT, Eigen::Dynamic, 1> const& center,
    RNGType& rng,
    unsigned int burnin = 0,
    SimplexBallWalkType walk_type = SimplexBallWalkType::GCW,
    double regcw_tau = 1.0,
    unsigned int max_reflections = 0)
{
    typedef typename Point::FT NT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;

    std::size_t const number_of_components = starting_points.size();

    if (number_of_components == 0)
    {
        throw std::invalid_argument(
            "run_simplex_ball_multichain: no starting points");
    }

    if (raw_component_weights.size() != number_of_components)
    {
        throw std::invalid_argument(
            "run_simplex_ball_multichain: "
            "number of weights must equal number of starting points");
    }

    unsigned int const dimension =
        static_cast<unsigned int>(A.cols());

    if (center.rows() != static_cast<int>(dimension))
    {
        throw std::invalid_argument(
            "run_simplex_ball_multichain: invalid center dimension");
    }

    for (VT const& starting_point : starting_points)
    {
        if (starting_point.rows() != static_cast<int>(dimension))
        {
            throw std::invalid_argument(
                "run_simplex_ball_multichain: "
                "invalid starting-point dimension");
        }
    }

    std::vector<std::size_t> assignments =
        sample_simplex_ball_component_assignments(
            raw_component_weights,
            number_of_samples,
            rng);

    std::vector<std::vector<unsigned int>> output_columns(
        number_of_components);

    for (unsigned int column = 0;
         column < number_of_samples;
         ++column)
    {
        output_columns[assignments[column]].push_back(column);
    }

    MT final_samples(dimension, number_of_samples);

    // One future and one independent RNG stream per active component.
    std::vector<std::future<MT>> component_futures(
        number_of_components);

    std::vector<bool> active_components(
        number_of_components,
        false);

    for (std::size_t component = 0;
        component < number_of_components;
        ++component)
    {
        unsigned int const component_sample_count =
            static_cast<unsigned int>(
                output_columns[component].size());

        if (component_sample_count == 0)
        {
            continue;
        }

        NT const seed_draw = rng.sample_urdist();

        unsigned int component_seed =
            static_cast<unsigned int>(
                seed_draw *
                static_cast<NT>(
                    std::numeric_limits<unsigned int>::max()));

        component_seed ^=
            0x9e3779b9u +
            static_cast<unsigned int>(component);

        RNGType component_rng(dimension);
        component_rng.set_seed(component_seed);

        active_components[component] = true;

        component_futures[component] =
            std::async(
                std::launch::async,
                [&, component, component_sample_count, component_rng]()
                    mutable -> MT
                {
                    return run_simplex_ball_chain<
                        Point,
                        MT,
                        RNGType>(
                            A,
                            b,
                            V,
                            starting_points[component],
                            component_sample_count,
                            walk_length,
                            center,
                            component_rng,
                            burnin,
                            walk_type,
                            regcw_tau,
                            max_reflections);
                });
    }

    // Merge on the calling thread. This avoids concurrent writes to the
    // final sample matrix and propagates worker exceptions through get().
    for (std::size_t component = 0;
        component < number_of_components;
        ++component)
    {
        if (!active_components[component])
        {
            continue;
        }

        MT component_samples =
            component_futures[component].get();

        unsigned int const component_sample_count =
            static_cast<unsigned int>(
                output_columns[component].size());

        for (unsigned int local_column = 0;
            local_column < component_sample_count;
            ++local_column)
        {
            unsigned int const output_column =
                output_columns[component][local_column];

            final_samples.col(output_column) =
                component_samples.col(local_column);
        }
    }

    return final_samples;
}

#endif
