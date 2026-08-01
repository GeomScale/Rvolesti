// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

using namespace Rcpp;
using Eigen::SparseMatrix;
using Eigen::VectorXd;

// [[Rcpp::export]]
SEXP create_sparse_hpolytope(const SparseMatrix<double> &A,
                             const VectorXd &b) {

    // Get dimensions
    int m = A.rows();
    int n = A.cols();

    // Proof that sparse matrix arrived correctly
    return List::create(
        Named("rows") = m,
        Named("cols") = n,
        Named("nonzeros") = A.nonZeros(),
        Named("b_length") = b.size()
    );
}