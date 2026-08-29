// [[Rcpp::depends(RcppEigen, BH)]]

// VolEsti (volume computation and sampling library)

// Copyright (c) 2026 Maria Zaza

// Licensed under GNU LGPL.3, see LICENSE file

#include <Rcpp.h>
#include <RcppEigen.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

#include <boost/random.hpp>

#include "generators/boost_random_number_generator.hpp"

namespace
{

typedef double NT;
typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;
typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;
typedef BoostRandomNumberGenerator<boost::mt19937, NT> RNGType;

bool segment_stays_outside_unit_ball(
    VT const& point,
    VT const& vertex,
    VT const& center,
    NT tolerance)
{
    VT const shifted_point = point - center;
    VT const direction = vertex - point;
    NT const squared_direction = direction.squaredNorm();

    if (squared_direction <= tolerance)
    {
        return shifted_point.squaredNorm() >= NT(1) - tolerance;
    }

    NT parameter =
        -shifted_point.dot(direction) / squared_direction;

    parameter = std::max(NT(0), std::min(NT(1), parameter));

    VT const nearest = shifted_point + parameter * direction;

    return nearest.squaredNorm() >= NT(1) - tolerance;
}

std::size_t nearest_starting_point(
    VT const& point,
    MT const& starting_points)
{
    std::size_t nearest = 0;
    NT best_distance = std::numeric_limits<NT>::infinity();

    for (Eigen::Index component = 0;
         component < starting_points.cols();
         ++component)
    {
        NT const distance =
            (point - starting_points.col(component)).squaredNorm();

        if (distance < best_distance)
        {
            best_distance = distance;
            nearest = static_cast<std::size_t>(component);
        }
    }

    return nearest;
}

} // namespace

//' Estimate relative simplex-sphere component surface areas
//'
//' This is an internal helper for `sample_simplex_ball()`. It uses uniform
//' directions on the complete sphere, keeps directions satisfying the
//' simplex inequalities, and classifies accepted points by the component of
//' the simplex 1-skeleton visible without entering the open unit ball. The
//' nearest validated starting point resolves only numerical or tangential
//' ambiguities.
//'
//' @param A,b Finite simplex half-space representation `A * x <= b`.
//' @param vertices Simplex vertices stored column-wise.
//' @param components List of one-based vertex-index components.
//' @param starting_points Matrix with one validated start per component.
//' @param center Centre of the unit sphere.
//' @param accepted_draws Requested number of feasible sphere draws.
//' @param max_attempts Maximum number of full-sphere proposals.
//' @param seed Optional non-negative RNG seed.
//' @param tolerance Positive numerical tolerance.
//' @return A list containing normalized component weights, Monte Carlo
//'   standard errors and counts, proposal metadata, ambiguous-assignment
//'   count, and any attempt-cap warning.
//'
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_simplex_ball_component_weights(
    Rcpp::NumericMatrix A,
    Rcpp::NumericVector b,
    Rcpp::NumericMatrix vertices,
    Rcpp::List components,
    Rcpp::NumericMatrix starting_points,
    Rcpp::NumericVector center,
    unsigned int accepted_draws,
    unsigned int max_attempts,
    Rcpp::Nullable<double> seed = R_NilValue,
    double tolerance = 1e-10)
{
    MT const matrix_a = Rcpp::as<MT>(A);
    VT const vector_b = Rcpp::as<VT>(b);
    MT const matrix_vertices = Rcpp::as<MT>(vertices);
    MT const matrix_starts = Rcpp::as<MT>(starting_points);
    VT const vector_center = Rcpp::as<VT>(center);

    Eigen::Index const dimension = matrix_a.cols();
    std::size_t const number_of_components =
        static_cast<std::size_t>(components.size());

    if (dimension < 2 ||
        matrix_a.rows() <= 0 ||
        vector_b.rows() != matrix_a.rows() ||
        matrix_vertices.rows() != dimension ||
        matrix_vertices.cols() <= 0 ||
        matrix_starts.rows() != dimension ||
        matrix_starts.cols() !=
            static_cast<Eigen::Index>(number_of_components) ||
        vector_center.rows() != dimension ||
        number_of_components == 0)
    {
        throw std::invalid_argument(
            "component-weight estimator received inconsistent dimensions");
    }

    if (!matrix_a.allFinite() ||
        !vector_b.allFinite() ||
        !matrix_vertices.allFinite() ||
        !matrix_starts.allFinite() ||
        !vector_center.allFinite())
    {
        throw std::invalid_argument(
            "component-weight estimator requires finite inputs");
    }

    if (accepted_draws == 0 || max_attempts < accepted_draws)
    {
        throw std::invalid_argument(
            "accepted draws must be positive and not exceed max attempts");
    }

    if (!std::isfinite(tolerance) || tolerance <= NT(0))
    {
        throw std::invalid_argument(
            "component-weight tolerance must be positive and finite");
    }

    std::vector<std::vector<Eigen::Index>> component_vertices;
    component_vertices.reserve(number_of_components);

    for (std::size_t component = 0;
         component < number_of_components;
         ++component)
    {
        Rcpp::IntegerVector const indices = components[component];

        if (indices.size() == 0)
        {
            throw std::invalid_argument(
                "every component must contain at least one vertex");
        }

        std::vector<Eigen::Index> converted;
        converted.reserve(static_cast<std::size_t>(indices.size()));

        for (int index : indices)
        {
            if (index == NA_INTEGER ||
                index < 1 ||
                index > matrix_vertices.cols())
            {
                throw std::invalid_argument(
                    "component vertex indices must be valid one-based indices");
            }

            converted.push_back(static_cast<Eigen::Index>(index - 1));
        }

        component_vertices.push_back(std::move(converted));
    }

    RNGType rng(static_cast<unsigned int>(dimension));

    if (seed.isNotNull())
    {
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

    std::vector<unsigned int> counts(number_of_components, 0);
    unsigned int accepted = 0;
    unsigned int attempts = 0;
    unsigned int fallback_assignments = 0;

    while (accepted < accepted_draws && attempts < max_attempts)
    {
        ++attempts;

        VT direction(dimension);
        NT norm = NT(0);

        do
        {
            for (Eigen::Index coordinate = 0;
                 coordinate < dimension;
                 ++coordinate)
            {
                direction(coordinate) = rng.sample_ndist();
            }

            norm = direction.norm();
        }
        while (!std::isfinite(norm) ||
               norm <= std::numeric_limits<NT>::epsilon());

        direction /= norm;
        VT const point = vector_center + direction;

        if ((matrix_a * point - vector_b).maxCoeff() > tolerance)
        {
            continue;
        }

        ++accepted;

        std::vector<std::size_t> visible_components;

        for (std::size_t component = 0;
             component < number_of_components;
             ++component)
        {
            bool visible = false;

            for (Eigen::Index vertex : component_vertices[component])
            {
                if (segment_stays_outside_unit_ball(
                        point,
                        matrix_vertices.col(vertex),
                        vector_center,
                        tolerance))
                {
                    visible = true;
                    break;
                }
            }

            if (visible)
            {
                visible_components.push_back(component);
            }
        }

        std::size_t assigned_component;

        if (visible_components.size() == 1)
        {
            assigned_component = visible_components.front();
        }
        else
        {
            assigned_component =
                nearest_starting_point(point, matrix_starts);
            ++fallback_assignments;
        }

        ++counts[assigned_component];
    }

    if (accepted == 0)
    {
        throw std::runtime_error(
            "component-weight estimator found no feasible sphere points");
    }

    Rcpp::NumericVector weights(number_of_components);
    Rcpp::NumericVector standard_errors(number_of_components);
    Rcpp::IntegerVector returned_counts(number_of_components);

    for (std::size_t component = 0;
         component < number_of_components;
         ++component)
    {
        NT const weight =
            static_cast<NT>(counts[component]) /
            static_cast<NT>(accepted);

        weights[component] = weight;
        standard_errors[component] =
            std::sqrt(
                weight * (NT(1) - weight) /
                static_cast<NT>(accepted));
        returned_counts[component] =
            static_cast<int>(counts[component]);
    }

    std::string warning;

    if (accepted < accepted_draws)
    {
        warning =
            "maximum attempts reached before the requested number of "
            "feasible sphere draws";
    }

    return Rcpp::List::create(
        Rcpp::Named("weights") = weights,
        Rcpp::Named("standard_errors") = standard_errors,
        Rcpp::Named("counts") = returned_counts,
        Rcpp::Named("accepted") = accepted,
        Rcpp::Named("attempts") = attempts,
        Rcpp::Named("fallback_assignments") = fallback_assignments,
        Rcpp::Named("warning") = warning);
}
