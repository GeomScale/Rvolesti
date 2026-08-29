// [[Rcpp::depends(RcppEigen, BH)]]

// VolEsti (volume computation and sampling library)

// Copyright (c) 2026 Maria Zaza

// Licensed under GNU LGPL.3, see LICENSE file

#include <Rcpp.h>
#include <RcppEigen.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <future>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/random.hpp>

#include "cartesian_geom/cartesian_kernel.h"
#include "generators/boost_random_number_generator.hpp"
#include "random_walks/simplex_ball_convergence.hpp"

namespace
{

typedef double NT;
typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;
typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;
typedef Cartesian<NT> Kernel;
typedef Kernel::Point Point;
typedef BoostRandomNumberGenerator<boost::mt19937, NT> RNGType;
typedef SimplexBallConvergenceResult<NT, VT, MT> ConvergenceResult;

Rcpp::List wrap_vector_list(std::vector<VT> const& vectors)
{
    Rcpp::List output(vectors.size());

    for (std::size_t index = 0; index < vectors.size(); ++index)
    {
        output[index] = Rcpp::wrap(vectors[index]);
    }

    return output;
}

Rcpp::List wrap_chain_list(std::vector<MT> const& chains)
{
    Rcpp::List output(chains.size());

    for (std::size_t chain = 0; chain < chains.size(); ++chain)
    {
        output[chain] = Rcpp::wrap(chains[chain]);
    }

    return output;
}

template <typename Diagnostics>
Rcpp::List wrap_diagnostics(
    Diagnostics const& diagnostics,
    unsigned int draws_per_chain)
{
    return Rcpp::List::create(
        Rcpp::Named("draws_per_chain") = draws_per_chain,
        Rcpp::Named("marginal_psrf") =
            Rcpp::wrap(diagnostics.marginal_psrf),
        Rcpp::Named("multivariate_psrf") =
            diagnostics.multivariate_psrf,
        Rcpp::Named("ess") =
            Rcpp::wrap(diagnostics.effective_sample_size),
        Rcpp::Named("per_chain_ess") =
            wrap_vector_list(diagnostics.per_chain_ess),
        Rcpp::Named("converged") = diagnostics.converged);
}

SimplexBallWalkType parse_walk_type(std::string const& walk)
{
    if (walk == "GCW")
    {
        return SimplexBallWalkType::GCW;
    }

    if (walk == "ReGCW")
    {
        return SimplexBallWalkType::ReGCW;
    }

    throw std::invalid_argument(
        "simplex-ball walk must be either 'GCW' or 'ReGCW'");
}

void apply_optional_seed(
    RNGType& rng,
    Rcpp::Nullable<double> const& seed)
{
    if (seed.isNull())
    {
        return;
    }

    NT const seed_value = Rcpp::as<NT>(seed);

    if (!std::isfinite(seed_value) ||
        seed_value < NT(0) ||
        seed_value >
            static_cast<NT>(
                std::numeric_limits<unsigned int>::max()))
    {
        throw std::invalid_argument(
            "seed must be finite and lie in the unsigned integer range");
    }

    rng.set_seed(static_cast<unsigned int>(seed_value));
}

Rcpp::List wrap_convergence_result(
    ConvergenceResult const& result,
    unsigned int check_interval)
{
    Rcpp::List history(result.diagnostic_history.size());

    for (std::size_t check = 0;
         check < result.diagnostic_history.size();
         ++check)
    {
        std::size_t const cumulative =
            (check + 1) * static_cast<std::size_t>(check_interval);

        unsigned int const draws_at_check =
            static_cast<unsigned int>(
                std::min<std::size_t>(
                    cumulative,
                    result.draws_per_chain));

        history[check] = wrap_diagnostics(
            result.diagnostic_history[check],
            draws_at_check);
    }

    return Rcpp::List::create(
        Rcpp::Named("samples") = Rcpp::wrap(result.samples),
        Rcpp::Named("chains") = wrap_chain_list(result.chains),
        Rcpp::Named("marginal_psrf") =
            Rcpp::wrap(result.diagnostics.marginal_psrf),
        Rcpp::Named("multivariate_psrf") =
            result.diagnostics.multivariate_psrf,
        Rcpp::Named("ess") =
            Rcpp::wrap(result.diagnostics.effective_sample_size),
        Rcpp::Named("per_chain_ess") =
            wrap_vector_list(result.diagnostics.per_chain_ess),
        Rcpp::Named("history") = history,
        Rcpp::Named("checks") = result.diagnostic_checks,
        Rcpp::Named("draws_per_chain") = result.draws_per_chain,
        Rcpp::Named("converged") = result.converged,
        Rcpp::Named("warning") = result.warning);
}

std::vector<NT> normalize_component_weights(
    Rcpp::NumericVector const& component_weights,
    std::size_t number_of_components)
{
    if (component_weights.size() !=
        static_cast<R_xlen_t>(number_of_components))
    {
        throw std::invalid_argument(
            "number of component weights must equal number of starts");
    }

    std::vector<NT> weights(number_of_components);
    NT total = NT(0);

    for (std::size_t component = 0;
         component < number_of_components;
         ++component)
    {
        NT const weight = component_weights[component];

        if (!std::isfinite(weight) || weight < NT(0))
        {
            throw std::invalid_argument(
                "component weights must be finite and non-negative");
        }

        weights[component] = weight;
        total += weight;
    }

    if (!std::isfinite(total) || total <= NT(0))
    {
        throw std::invalid_argument(
            "component weights must have a positive finite sum");
    }

    for (NT& weight : weights)
    {
        weight /= total;
    }

    return weights;
}

std::vector<unsigned int> allocate_component_samples(
    unsigned int number_of_samples,
    std::vector<NT> const& weights)
{
    std::vector<unsigned int> counts(weights.size(), 0);
    std::vector<NT> remainders(weights.size(), NT(0));
    unsigned long long allocated = 0;

    for (std::size_t component = 0;
         component < weights.size();
         ++component)
    {
        NT const quota =
            static_cast<NT>(number_of_samples) * weights[component];
        unsigned int const count =
            static_cast<unsigned int>(std::floor(quota));

        counts[component] = count;
        remainders[component] = quota - static_cast<NT>(count);
        allocated += count;
    }

    if (allocated > number_of_samples)
    {
        throw std::logic_error(
            "largest-remainder component allocation exceeded the sample total");
    }

    unsigned int const remaining =
        number_of_samples - static_cast<unsigned int>(allocated);
    std::vector<std::size_t> order(weights.size());
    std::iota(order.begin(), order.end(), std::size_t(0));
    std::stable_sort(
        order.begin(),
        order.end(),
        [&](std::size_t left, std::size_t right)
        {
            return remainders[left] > remainders[right];
        });

    if (remaining > order.size())
    {
        throw std::logic_error(
            "largest-remainder component allocation failed");
    }

    for (unsigned int index = 0; index < remaining; ++index)
    {
        ++counts[order[index]];
    }

    return counts;
}

std::vector<unsigned int> sample_columns_without_replacement(
    unsigned int available,
    unsigned int requested,
    RNGType& rng)
{
    if (requested > available)
    {
        throw std::runtime_error(
            "convergence run retained fewer points than allocated");
    }

    std::vector<unsigned int> indices(available);
    std::iota(indices.begin(), indices.end(), 0u);

    for (unsigned int position = 0;
         position < requested;
         ++position)
    {
        unsigned int const remaining = available - position;
        unsigned int offset = static_cast<unsigned int>(
            rng.sample_urdist() * static_cast<NT>(remaining));

        if (offset >= remaining)
        {
            offset = remaining - 1;
        }

        std::swap(indices[position], indices[position + offset]);
    }

    indices.resize(requested);
    return indices;
}

} // namespace

//' Sample one connected simplex-sphere component until convergence
//'
//' This is the native convergence-controlled sampler used internally by
//' `sample_simplex_ball()`.
//'
//' @param A,b Finite simplex half-space representation `A * x <= b`.
//' @param V Simplex vertices stored column-wise.
//' @param start,center Starting point and unit-sphere centre.
//' @param n_chains Number of independent chains; at least two.
//' @param check_interval Retained draws per chain between diagnostic checks.
//' @param max_iterations Maximum retained draws per chain.
//' @param walk_length Number of transitions between retained draws.
//' @param burnin Number of discarded transitions per chain.
//' @param psrf_target,ess_target Convergence thresholds.
//' @param walk Either `"GCW"` or `"ReGCW"`.
//' @param regcw_tau Reflective-walk trajectory scale.
//' @param max_reflections Reflective-walk cap; zero selects the native default.
//' @param seed Optional non-negative RNG seed.
//' @return A list containing retained samples, separate chains, PSRF and ESS
//'   diagnostics, diagnostic history, convergence status, and any cap warning.
//'
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_sample_simplex_ball_component(
    Rcpp::NumericMatrix A,
    Rcpp::NumericVector b,
    Rcpp::NumericMatrix V,
    Rcpp::NumericVector start,
    Rcpp::NumericVector center,
    unsigned int n_chains,
    unsigned int check_interval,
    unsigned int max_iterations,
    unsigned int walk_length,
    unsigned int burnin,
    double psrf_target,
    double ess_target,
    std::string walk,
    double regcw_tau,
    unsigned int max_reflections,
    Rcpp::Nullable<double> seed = R_NilValue)
{
    MT const matrix_a = Rcpp::as<MT>(A);
    VT const vector_b = Rcpp::as<VT>(b);
    MT const matrix_vertices = Rcpp::as<MT>(V);
    VT const vector_start = Rcpp::as<VT>(start);
    VT const vector_center = Rcpp::as<VT>(center);

    Eigen::Index const dimension = matrix_a.cols();

    if (dimension < 2 ||
        matrix_a.rows() <= 0 ||
        vector_b.rows() != matrix_a.rows() ||
        matrix_vertices.rows() != dimension ||
        matrix_vertices.cols() <= 0 ||
        vector_start.rows() != dimension ||
        vector_center.rows() != dimension)
    {
        throw std::invalid_argument(
            "simplex-ball sampler received inconsistent dimensions");
    }

    if (!matrix_a.allFinite() ||
        !vector_b.allFinite() ||
        !matrix_vertices.allFinite() ||
        !vector_start.allFinite() ||
        !vector_center.allFinite())
    {
        throw std::invalid_argument(
            "simplex-ball sampler requires finite inputs");
    }

    SimplexBallWalkType const walk_type = parse_walk_type(walk);

    RNGType rng(static_cast<int>(dimension));
    apply_optional_seed(rng, seed);

    ConvergenceResult const result =
        run_simplex_ball_chains_until_converged<Point, MT, RNGType>(
            matrix_a,
            vector_b,
            matrix_vertices,
            vector_start,
            n_chains,
            check_interval,
            max_iterations,
            walk_length,
            vector_center,
            rng,
            burnin,
            psrf_target,
            ess_target,
            walk_type,
            regcw_tau,
            max_reflections);

    return wrap_convergence_result(result, check_interval);
}

//' Sample every active simplex-sphere component in parallel
//'
//' This is the native multicomponent sampler used by the high-level
//' `sample_simplex_ball()` pipeline. Each active component receives its own
//' RNG stream and convergence-controlled set of chains. Components are run
//' concurrently; PSRF and ESS are still evaluated only among chains targeting
//' the same connected component.
//'
//' @param A,b Finite simplex half-space representation `A * x <= b`.
//' @param V Simplex vertices stored column-wise.
//' @param starting_points Matrix with one feasible starting point per column.
//' @param component_weights Non-negative component surface-area weights.
//' @param N Exact total number of samples to return. The native wrapper uses
//'   deterministic largest-remainder allocation across components.
//' @param center Centre of the unit sphere.
//' @param n_chains Number of independent diagnostic chains per component.
//' @param check_interval Base retained draws per chain between diagnostic
//'   checks. It is increased when necessary to provide each component pool.
//' @param max_iterations Maximum retained draws per chain.
//' @param walk_length Number of transitions between retained draws.
//' @param burnin Number of discarded transitions per chain.
//' @param psrf_target,ess_target Convergence thresholds.
//' @param walk Either `"GCW"` or `"ReGCW"`.
//' @param regcw_tau Reflective-walk trajectory scale.
//' @param max_reflections Reflective-walk cap; zero selects the native default.
//' @param seed Optional non-negative master RNG seed.
//' @return A list containing the exact combined sample matrix, one-based
//'   sample-component labels, normalized weights and counts, per-component
//'   convergence results, aggregate PSRF/ESS, convergence status and warnings.
//'
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_sample_simplex_ball(
    Rcpp::NumericMatrix A,
    Rcpp::NumericVector b,
    Rcpp::NumericMatrix V,
    Rcpp::NumericMatrix starting_points,
    Rcpp::NumericVector component_weights,
    Rcpp::NumericVector center,
    unsigned int N,
    unsigned int n_chains,
    unsigned int check_interval,
    unsigned int max_iterations,
    unsigned int walk_length,
    unsigned int burnin,
    double psrf_target,
    double ess_target,
    std::string walk,
    double regcw_tau,
    unsigned int max_reflections,
    Rcpp::Nullable<double> seed = R_NilValue)
{
    MT const matrix_a = Rcpp::as<MT>(A);
    VT const vector_b = Rcpp::as<VT>(b);
    MT const matrix_vertices = Rcpp::as<MT>(V);
    MT const matrix_starts = Rcpp::as<MT>(starting_points);
    VT const vector_center = Rcpp::as<VT>(center);

    Eigen::Index const dimension = matrix_a.cols();
    std::size_t const number_of_components =
        static_cast<std::size_t>(matrix_starts.cols());

    if (dimension < 2 ||
        matrix_a.rows() <= 0 ||
        vector_b.rows() != matrix_a.rows() ||
        matrix_vertices.rows() != dimension ||
        matrix_vertices.cols() <= 0 ||
        matrix_starts.rows() != dimension ||
        number_of_components == 0 ||
        vector_center.rows() != dimension)
    {
        throw std::invalid_argument(
            "parallel simplex-ball sampler received inconsistent dimensions");
    }

    if (!matrix_a.allFinite() ||
        !vector_b.allFinite() ||
        !matrix_vertices.allFinite() ||
        !matrix_starts.allFinite() ||
        !vector_center.allFinite())
    {
        throw std::invalid_argument(
            "parallel simplex-ball sampler requires finite inputs");
    }

    if (N == 0 || n_chains < 2)
    {
        throw std::invalid_argument(
            "parallel simplex-ball sampling needs samples and two chains");
    }

    std::vector<NT> const normalized_weights =
        normalize_component_weights(
            component_weights,
            number_of_components);
    std::vector<unsigned int> const counts =
        allocate_component_samples(N, normalized_weights);

    for (std::size_t component = 0;
         component < number_of_components;
         ++component)
    {
        if (static_cast<unsigned long long>(counts[component]) >
            static_cast<unsigned long long>(n_chains) *
                static_cast<unsigned long long>(max_iterations))
        {
            throw std::invalid_argument(
                "maximum iterations cannot provide the component allocation");
        }

    }

    SimplexBallWalkType const walk_type = parse_walk_type(walk);
    RNGType master_rng(static_cast<int>(dimension));
    apply_optional_seed(master_rng, seed);

    std::vector<std::future<ConvergenceResult>> futures(
        number_of_components);
    std::vector<unsigned int> component_intervals(
        number_of_components,
        0);
    std::vector<bool> active(number_of_components, false);

    for (std::size_t component = 0;
         component < number_of_components;
         ++component)
    {
        if (counts[component] == 0)
        {
            continue;
        }

        unsigned int const required_draws =
            static_cast<unsigned int>(
                (static_cast<unsigned long long>(counts[component]) +
                 static_cast<unsigned long long>(n_chains) - 1) /
                static_cast<unsigned long long>(n_chains));
        unsigned int const component_interval =
            std::max(
                std::max(check_interval, 4u),
                required_draws);

        NT const seed_draw = master_rng.sample_urdist();
        unsigned int component_seed =
            static_cast<unsigned int>(
                seed_draw *
                static_cast<NT>(
                    std::numeric_limits<unsigned int>::max()));
        component_seed ^=
            0x9e3779b9u + static_cast<unsigned int>(component);

        VT const component_start =
            matrix_starts.col(static_cast<Eigen::Index>(component));

        component_intervals[component] = component_interval;
        active[component] = true;

        futures[component] = std::async(
            std::launch::async,
            [&, component_start, component_interval, component_seed]()
                mutable -> ConvergenceResult
            {
                RNGType component_rng(static_cast<int>(dimension));
                component_rng.set_seed(component_seed);

                return run_simplex_ball_chains_until_converged<
                    Point,
                    MT,
                    RNGType>(
                        matrix_a,
                        vector_b,
                        matrix_vertices,
                        component_start,
                        n_chains,
                        component_interval,
                        max_iterations,
                        walk_length,
                        vector_center,
                        component_rng,
                        burnin,
                        psrf_target,
                        ess_target,
                        walk_type,
                        regcw_tau,
                        max_reflections);
            });
    }

    Rcpp::List component_results(number_of_components);
    Rcpp::LogicalVector active_components(number_of_components);
    Rcpp::IntegerVector returned_counts(number_of_components);
    MT grouped_samples(dimension, static_cast<Eigen::Index>(N));
    std::vector<unsigned int> grouped_components(N);
    unsigned int output_column = 0;
    NT overall_psrf = -std::numeric_limits<NT>::infinity();
    NT total_ess = NT(0);
    bool all_converged = true;
    std::vector<std::string> warnings;

    // Only this calling thread touches R objects. Worker exceptions are
    // propagated by future::get() before any partial result is returned.
    for (std::size_t component = 0;
         component < number_of_components;
         ++component)
    {
        active_components[component] = active[component];
        returned_counts[component] =
            static_cast<int>(counts[component]);

        if (!active[component])
        {
            component_results[component] = R_NilValue;
            continue;
        }

        ConvergenceResult const result = futures[component].get();

        std::vector<unsigned int> const selected =
            sample_columns_without_replacement(
                static_cast<unsigned int>(result.samples.cols()),
                counts[component],
                master_rng);

        for (unsigned int local_column = 0;
             local_column < counts[component];
             ++local_column)
        {
            grouped_samples.col(output_column) =
                result.samples.col(selected[local_column]);
            grouped_components[output_column] =
                static_cast<unsigned int>(component + 1);
            ++output_column;
        }

        NT component_psrf = result.diagnostics.multivariate_psrf;

        if (result.diagnostics.marginal_psrf.rows() > 0)
        {
            component_psrf = std::max(
                component_psrf,
                result.diagnostics.marginal_psrf.maxCoeff());
        }

        overall_psrf = std::max(overall_psrf, component_psrf);
        total_ess +=
            result.diagnostics.effective_sample_size.minCoeff();
        all_converged = all_converged && result.converged;

        if (!result.warning.empty())
        {
            warnings.push_back(
                "Component " + std::to_string(component + 1) +
                ": " + result.warning);
        }

        component_results[component] = wrap_convergence_result(
            result,
            component_intervals[component]);
    }

    if (output_column != N)
    {
        throw std::logic_error(
            "parallel simplex-ball sampler produced an invalid allocation");
    }

    std::vector<unsigned int> const shuffled_columns =
        sample_columns_without_replacement(N, N, master_rng);
    MT samples(dimension, static_cast<Eigen::Index>(N));
    Rcpp::IntegerVector sample_components(N);

    for (unsigned int column = 0; column < N; ++column)
    {
        unsigned int const source = shuffled_columns[column];
        samples.col(column) = grouped_samples.col(source);
        sample_components[column] =
            static_cast<int>(grouped_components[source]);
    }

    return Rcpp::List::create(
        Rcpp::Named("samples") = Rcpp::wrap(samples),
        Rcpp::Named("sample_components") = sample_components,
        Rcpp::Named("component_weights") =
            Rcpp::wrap(normalized_weights),
        Rcpp::Named("component_counts") = returned_counts,
        Rcpp::Named("component_results") = component_results,
        Rcpp::Named("active_components") = active_components,
        Rcpp::Named("psrf") = overall_psrf,
        Rcpp::Named("ess") = total_ess,
        Rcpp::Named("converged") = all_converged,
        Rcpp::Named("warning") = Rcpp::wrap(warnings),
        Rcpp::Named("parallel_components") = true);
}
