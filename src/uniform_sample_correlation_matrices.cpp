#include <Rcpp.h>
#include <RcppEigen.h>
#include <iostream>
#include <chrono>
#include <vector>
#include <boost/random.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include "cartesian_geom/cartesian_kernel.h"
#include "random_walks/random_walks.hpp"
#include "sampling/sample_correlation_matrices.hpp"
#include "matrix_operations/EigenvaluesProblems.h"

//' Uniformly sample correlation matrices
//'
//' @param n The dimension of the correlation matrix.
//' @param num_matrices The number of matrices to generate.
//' @param walk_length The length of the random walk.
//' @param nburns The number of burn-in steps for the random walk.
//' @param validate Optional. Whether to validate the sampled matrices. Default is false.
//'
//' @return A list of sampled correlation matrices.

// [[Rcpp::export]]
Rcpp::List uniform_sample_correlation_matrices(const unsigned int n, const unsigned int num_matrices = 1000,
                                               const unsigned int walk_length=1, const unsigned int nburns = 0,  const bool validate = false) {
    typedef double NT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;
    typedef BoostRandomNumberGenerator<boost::mt19937, NT> RNGType;
    typedef CorreMatrix<NT> PointMT;

    std::list<MT> randCorMatrices;

    uniform_correlation_sampling_MT<AcceleratedBilliardWalk, PointMT, RNGType>(n, randCorMatrices, walk_length, num_matrices, nburns);

    EigenvaluesProblems<NT, MT, Eigen::Matrix<NT, Eigen::Dynamic, 1>> solver;

     if (validate) {
         const NT tol = 1e-8;
         for (const auto& mat : randCorMatrices) {
             if (!solver.is_correlation_matrix(mat, tol)) {
                 throw Rcpp::exception("Invalid correlation matrix");
             }
         }
     } 

    Rcpp::List rcpp_sampled_matrices(randCorMatrices.size());
    int index = 0;
    for (const auto& mat : randCorMatrices) {
        rcpp_sampled_matrices[index++] = Rcpp::wrap(mat);
    }

    Rcpp::List result = Rcpp::List::create(
        Rcpp::Named("sampled_matrices") = rcpp_sampled_matrices
    );

    return result;
}
