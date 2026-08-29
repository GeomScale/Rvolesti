// [[Rcpp::depends(RcppEigen, BH)]]

// VolEsti (volume computation and sampling library)

// Copyright (c) 2026 Maria Zaza

// Licensed under GNU LGPL.3, see LICENSE file

#include <Rcpp.h>
#include <RcppEigen.h>

#include <cmath>
#include <stdexcept>
#include <vector>

#include "convex_bodies/simplexintersectball_components.h"

namespace
{

typedef double NT;
typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;
typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;

} // namespace

//' Find connected components of a simplex-sphere intersection
//'
//' This is an internal helper for `sample_simplex_ball()`. Vertex indices in
//' the returned components are one-based, as expected by R.
//'
//' @param vertices Simplex vertices stored column-wise.
//' @param center Centre of the ball.
//' @param radius Positive ball radius.
//' @param tol Positive numerical tolerance.
//' @return A list containing the one-based components, active-vertex mask,
//'   active vertex indices, and 1-skeleton adjacency matrix.
//'
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_simplex_ball_components(
    Rcpp::NumericMatrix vertices,
    Rcpp::NumericVector center,
    double radius = 1.0,
    double tol = 1e-10)
{
    MT const matrix_vertices = Rcpp::as<MT>(vertices);
    VT const vector_center = Rcpp::as<VT>(center);

    if (matrix_vertices.rows() <= 0 ||
        matrix_vertices.cols() <= 0 ||
        vector_center.rows() != matrix_vertices.rows())
    {
        throw std::invalid_argument(
            "simplex-ball components received inconsistent dimensions");
    }

    if (!matrix_vertices.allFinite() || !vector_center.allFinite())
    {
        throw std::invalid_argument(
            "simplex-ball components require finite vertices and center");
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

    std::vector<int> const active =
        active_vertices_outside_ball(
            matrix_vertices,
            vector_center,
            radius,
            tol);

    Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> const adjacency =
        build_simplex_ball_graph(
            matrix_vertices,
            vector_center,
            radius,
            tol);

    std::vector<std::vector<int>> const components =
        connected_components_from_graph(adjacency, active);

    Rcpp::List component_list(components.size());

    for (std::size_t component = 0;
         component < components.size();
         ++component)
    {
        Rcpp::IntegerVector indices(components[component].size());

        for (std::size_t index = 0;
             index < components[component].size();
             ++index)
        {
            indices[index] = components[component][index] + 1;
        }

        component_list[component] = indices;
    }

    Rcpp::LogicalVector active_mask(active.size());
    std::vector<int> active_indices;
    active_indices.reserve(active.size());

    for (std::size_t vertex = 0; vertex < active.size(); ++vertex)
    {
        bool const is_active = active[vertex] != 0;
        active_mask[vertex] = is_active;

        if (is_active)
        {
            active_indices.push_back(static_cast<int>(vertex) + 1);
        }
    }

    return Rcpp::List::create(
        Rcpp::Named("components") = component_list,
        Rcpp::Named("active") = active_mask,
        Rcpp::Named("active_vertices") = Rcpp::wrap(active_indices),
        Rcpp::Named("adjacency") = Rcpp::wrap(adjacency));
}
