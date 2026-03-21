// [[Rcpp::depends(BH)]]

// VolEsti (volume computation and sampling library)
// Licensed under GNU LGPL.3, see LICENCE file

// Count linear extensions of a partial order via volume estimation
// of the corresponding order polytope.

#include <Rcpp.h>
#include <RcppEigen.h>
#include <boost/random.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>

#include "cartesian_geom/cartesian_kernel.h"
#include "convex_bodies/orderpolytope.h"
#include "misc/poset.h"
#include "random_walks/random_walks.hpp"
#include "volume/volume_sequence_of_balls.hpp"
#include "volume/volume_cooling_gaussians.hpp"
#include "volume/volume_cooling_balls.hpp"


//' Count the number of linear extensions of a partial order
//'
//' Uses volume approximation of the corresponding order polytope.
//' The number of linear extensions equals n! * Vol(OrderPolytope).
//'
//' @param n Integer. Number of elements in the poset (elements are 0 to n-1).
//' @param relations A matrix with 2 columns. Each row (i, j) represents the relation a_i <= a_j.
//' @param algorithm String. Volume algorithm to use: "CB" (Cooling Balls, default), "CG" (Cooling Gaussians), or "SOB" (Sequence of Balls).
//' @param error Numeric. Upper bound for volume approximation error (default: 0.1).
//' @param walk_length Integer. Walk length for the random walk (default: 1).
//'
//' @return A list with:
//' \describe{
//'   \item{linear_extensions}{Approximate number of linear extensions}
//'   \item{volume}{Approximate volume of the order polytope}
//'   \item{dimension}{Number of elements in the poset}
//' }
//'
//' @examples
//' # Chain: a0 <= a1 <= a2 <= a3 (exactly 1 linear extension)
//' result <- count_linear_extensions(4, rbind(c(0,1), c(1,2), c(2,3)))
//'
//' # Antichain: no relations (exactly 4! = 24 linear extensions)
//' result <- count_linear_extensions(4, matrix(ncol=2, nrow=0))
//'
//' @export
// [[Rcpp::export]]
Rcpp::List count_linear_extensions(int n,
                                    Rcpp::Nullable<Rcpp::IntegerMatrix> relations = R_NilValue,
                                    Rcpp::Nullable<std::string> algorithm = R_NilValue,
                                    double error = 0.1,
                                    int walk_length = 1)
{
    typedef double NT;
    typedef Cartesian<NT> Kernel;
    typedef typename Kernel::Point Point;
    typedef BoostRandomNumberGenerator<boost::mt19937, NT> RNGType;
    typedef OrderPolytope<Point> OrderPoly;
    typedef typename OrderPoly::VT VT;
    typedef typename OrderPoly::MT MT;

    // Build the poset
    typedef typename Poset::RT RT;
    typedef typename Poset::RV RV;

    RV order_relations;
    if (relations.isNotNull()) {
        Rcpp::IntegerMatrix rel_mat(relations);
        for (int i = 0; i < rel_mat.nrow(); ++i) {
            order_relations.push_back(RT(rel_mat(i, 0), rel_mat(i, 1)));
        }
    }

    Poset poset(n, order_relations);
    OrderPoly OP(poset);

    unsigned int d = OP.dimension();
    RNGType rng(d);

    // Determine algorithm
    std::string algo = "CB";
    if (algorithm.isNotNull()) {
        algo = Rcpp::as<std::string>(algorithm);
    }

    // Compute volume
    NT volume;
    if (algo == "SOB" || algo == "sob") {
        unsigned int wl = (walk_length == 1) ? (10 + d / 10) : walk_length;
        volume = volume_sequence_of_balls<CDHRWalk, RNGType>(OP, rng, error, wl);
    } else if (algo == "CG" || algo == "cg") {
        volume = volume_cooling_gaussians<GaussianCDHRWalk, RNGType>(OP, rng, error, walk_length);
    } else {
        volume = volume_cooling_balls<CDHRWalk, RNGType>(OP, rng, error, walk_length).second;
    }

    // Compute n! * volume
    NT le_count = volume;
    for (int i = d; i > 1; --i) {
        le_count *= NT(i);
    }

    return Rcpp::List::create(
        Rcpp::Named("linear_extensions") = le_count,
        Rcpp::Named("volume") = volume,
        Rcpp::Named("dimension") = d
    );
}
