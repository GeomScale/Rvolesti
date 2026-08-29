// [[Rcpp::depends(RcppEigen, BH)]]

// VolEsti (volume computation and sampling library)

// Copyright (c) 2026 Maria Zaza

// Licensed under GNU LGPL.3, see LICENSE file

#include <Rcpp.h>
#include <RcppEigen.h>

#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>

#include "convex_bodies/simplexintersectball_components.h"

namespace
{

typedef double NT;
typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;
typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;

void validate_interior_point(
    MT const& A,
    VT const& b,
    VT const& point,
    VT const& center,
    NT radius,
    NT tol)
{
    if (!point.allFinite())
    {
        throw std::invalid_argument(
            "simplex-ball interior point must be finite");
    }

    if (!point_is_inside_ball(point, center, radius, tol))
    {
        throw std::invalid_argument(
            "simplex-ball interior point must lie strictly inside the ball");
    }

    if (!point_satisfies_halfspaces(A, b, point, tol))
    {
        throw std::invalid_argument(
            "simplex-ball interior point violates the simplex inequalities");
    }
}

} // namespace

//' Find a starting point for one simplex-sphere component
//'
//' This is an internal helper for `sample_simplex_ball()`. `component` uses
//' one-based R vertex indices; the native component routine uses zero-based
//' indices internally.
//'
//' @param vertices Simplex vertices stored column-wise.
//' @param component One-based vertex indices for one connected component.
//' @param A,b Finite simplex half-space representation `A * x <= b`.
//' @param center Centre of the ball.
//' @param interior_point Optional strict interior point. When omitted, a
//'   Chebyshev-centre approximation is computed natively.
//' @param radius Positive ball radius.
//' @param tol Positive numerical tolerance.
//' @return A feasible point on the requested spherical component.
//'
//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericVector cpp_simplex_ball_start(
    Rcpp::NumericMatrix vertices,
    Rcpp::IntegerVector component,
    Rcpp::NumericMatrix A,
    Rcpp::NumericVector b,
    Rcpp::NumericVector center,
    Rcpp::Nullable<Rcpp::NumericVector> interior_point = R_NilValue,
    double radius = 1.0,
    double tol = 1e-10)
{
    MT const matrix_vertices = Rcpp::as<MT>(vertices);
    MT const matrix_a = Rcpp::as<MT>(A);
    VT const vector_b = Rcpp::as<VT>(b);
    VT const vector_center = Rcpp::as<VT>(center);

    Eigen::Index const dimension = matrix_vertices.rows();

    if (dimension <= 0 ||
        matrix_vertices.cols() <= 0 ||
        matrix_a.rows() <= 0 ||
        matrix_a.cols() != dimension ||
        vector_b.rows() != matrix_a.rows() ||
        vector_center.rows() != dimension)
    {
        throw std::invalid_argument(
            "simplex-ball starting-point search received inconsistent dimensions");
    }

    if (!matrix_vertices.allFinite() ||
        !matrix_a.allFinite() ||
        !vector_b.allFinite() ||
        !vector_center.allFinite())
    {
        throw std::invalid_argument(
            "simplex-ball starting-point search requires finite inputs");
    }

    if (component.size() == 0)
    {
        throw std::invalid_argument(
            "simplex-ball component must contain at least one vertex");
    }

    if (!std::isfinite(radius) || radius <= NT(0))
    {
        throw std::invalid_argument(
            "simplex-ball radius must be positive and finite");
    }

    if (!std::isfinite(tol) || tol <= NT(0))
    {
        throw std::invalid_argument(
            "simplex-ball tolerance must be positive and finite");
    }

    std::vector<int> native_component;
    native_component.reserve(static_cast<std::size_t>(component.size()));

    std::vector<bool> seen(
        static_cast<std::size_t>(matrix_vertices.cols()),
        false);

    for (int index : component)
    {
        if (index == NA_INTEGER ||
            index < 1 ||
            index > matrix_vertices.cols())
        {
            throw std::invalid_argument(
                "component vertex indices must be valid one-based indices");
        }

        std::size_t const zero_based =
            static_cast<std::size_t>(index - 1);

        if (seen[zero_based])
        {
            throw std::invalid_argument(
                "component vertex indices must not contain duplicates");
        }

        seen[zero_based] = true;
        native_component.push_back(index - 1);
    }

    VT interior(dimension);

    if (interior_point.isNotNull())
    {
        Rcpp::NumericVector const supplied(interior_point);

        if (supplied.size() != dimension)
        {
            throw std::invalid_argument(
                "simplex-ball interior point has the wrong dimension");
        }

        interior = Rcpp::as<VT>(supplied);
    }
    else
    {
        std::pair<bool, VT> const chebyshev_result =
            chebyshev_center_intersect_ball(
                matrix_a,
                vector_b,
                vector_center,
                radius,
                50,
                tol);

        if (!chebyshev_result.first)
        {
            Rcpp::stop(
                "unable to compute an interior point for the simplex-ball intersection");
        }

        interior = chebyshev_result.second;
    }

    validate_interior_point(
        matrix_a,
        vector_b,
        interior,
        vector_center,
        radius,
        tol);

    std::pair<bool, VT> const starting_result =
        find_starting_point_for_component(
            matrix_vertices,
            native_component,
            matrix_a,
            vector_b,
            interior,
            vector_center,
            radius,
            tol);

    if (!starting_result.first)
    {
        Rcpp::stop(
            "unable to find a feasible starting point for the requested simplex-ball component");
    }

    return Rcpp::wrap(starting_result.second);
}
