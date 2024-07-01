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

//' Uniformly sample correlation matrices
//'
//' @param n The dimension of the correlation matrix.
//' @param num_points The number of sample points to generate.
//' @param validate Optional. Whether to validate the sampled matrices. Default is false.
//'
//' @return A list of sampled correlation matrices.

// [[Rcpp::export]]
Rcpp::List uniform_sampling_correlation_MT(const unsigned int n, const unsigned int num_points = 1000, const bool validate = false) {
    typedef double NT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;
    typedef BoostRandomNumberGenerator<boost::mt19937, NT, 3> RNGType;
    typedef CorreMatrix<NT> PointMT;

    std::cout << num_points << " correlation matrices of size " << n << " with matrix PointType" << std::endl;
    std::chrono::steady_clock::time_point start, end_time;
    double time;
    std::vector<PointMT> randPoints;
    unsigned int walkL = 1;

    start = std::chrono::steady_clock::now();

    uniform_correlation_sampling_MT<AcceleratedBilliardWalk, PointMT, RNGType>(n, randPoints, walkL, num_points, 0);

    end_time = std::chrono::steady_clock::now();
    time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start).count();
    std::cout << "Elapsed time : " << time << " (ms)" << std::endl;

    int valid_points = 0;
    const double tol = 1e-8;

    std::vector<MT> sampled_points; 

    for (const auto& points : randPoints) {
        sampled_points.push_back(points.mat); // Store the sampled point

        if (validate) {
            bool flag = true;

            // Check if all the diagonal elements are 1
            for (int i = 0; i < points.mat.rows(); i++) {
                if (std::abs(points.mat(i, i) - 1.0) > tol) {
                    flag = false;
                    break;
                }
            }

            // Check if the matrix is positive definite
            Eigen::SelfAdjointEigenSolver<MT> eigen_solver(points.mat);

            if (eigen_solver.info() != Eigen::Success) flag = false;

            if (eigen_solver.eigenvalues().minCoeff() > -tol) flag = true;
            else flag = false;

            if (flag == true) valid_points++;
        }
    }

    if (validate) {
        std::cout << "Number of valid points = " << valid_points << std::endl;
    }

    Rcpp::List rcpp_sampled_points(sampled_points.size());
    for (size_t i = 0; i < sampled_points.size(); ++i) {
        rcpp_sampled_points[i] = Rcpp::wrap(sampled_points[i]);
    }

    Rcpp::List result = Rcpp::List::create(
        Rcpp::Named("elapsed_time") = time,
        Rcpp::Named("sampled_points") = rcpp_sampled_points
    );

    return result;
}
