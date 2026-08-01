#' Sample points from a convex polytope
#'
#' Sample perfect uniformly distributed points from well known convex bodies
#' or run geometric random walks to approximate arbitrary distributions.
#'
#' The \eqn{d}-dimensional unit simplex is the set of points \eqn{\vec{x} \in \R^d},
#' s.t.: \eqn{\sum_i x_i \leq 1}, \eqn{x_i \geq 0}. The \eqn{d}-dimensional canonical
#' simplex is the set of points \eqn{\vec{x} \in \R^d}, s.t.: \eqn{\sum_i x_i = 1},
#' \eqn{x_i \geq 0}.
#'
#' @param P A convex body object. One of: Hpolytope, Vpolytope, Zonotope,
#'   VpolytopeIntersection, or Spectrahedron.
#' @param n The number of points that the function is going to sample.
#' @param random_walk Optional. A list that declares the random walk method and
#'   related parameters as follows:
#'   \describe{
#'     \item{\code{walk}}{A string: (a) \code{'CDHR'} for coordinate direction hit-and-run,
#'       b) \code{'RDHR'} for random direction hit-and-run, c) \code{'BaW'} for the ball walk,
#'       d) \code{'BiW'} for the billiard walk, e) \code{'BiRCDHR'} for the bicriteria random walk,
#'       f) \code{'Dikin'} for dikin walk, g) \code{'Vaidya'} for vaidya walk,
#'       h) \code{'John'} for john walk, i) \code{'BRDHR'} for boundary hit-and-run,
#'       j) \code{'Hamiltonian'} for Hamiltonian walk, k) \code{'accelerated_billiard'} for
#'       accelerated billiard walk, l) \code{'BilliardRef'} for billiard walk with reflections,
#'       m) \code{'HRlu'} for logconcave settings with H-polytope,
#'       n) \code{'HRnr'} for logconcave densities with V-polytope, o) \code{'HESSIANRWR'}
#'       for logconcave densities with H-polytope and sparse constrained problems.}
#'     \item{\code{walk_length}}{The number of the steps per generated point for the random walk.
#'       The default value is \eqn{1}.}
#'     \item{\code{nburns}}{The number of points to burn before start sampling. The default
#'       value is \eqn{1}.}
#'     \item{\code{starting_point}}{A \eqn{d}-dimensional numerical vector that declares a
#'       starting point in the interior of the polytope for the random walk. The default choice
#'       is the center of the ball as computed by \code{inner_ball()}.}
#'     \item{\code{BaW_rad}}{The radius for the ball walk.}
#'     \item{\code{L}}{The maximum length of the billiard trajectory or the radius for the step
#'       of Dikin, Vaidya or John walk.}
#'     \item{\code{solver}}{Specify ODE solver for logconcave sampling. Options are i) leapfrog,
#'       ii) euler iii) runge-kutta iv) richardson}
#'     \item{\code{step_size}}{Optionally chosen step size for logconcave sampling. Defaults to
#'       a theoretical value if not provided.}
#'   }
#' @param distribution Optional. A list that declares the target density and related
#'   parameters as follows:
#'   \describe{
#'     \item{\code{density}}{A string: (a) \code{'uniform'} for the uniform distribution,
#'       (b) \code{'gaussian'} for the multidimensional spherical distribution,
#'       (c) \code{'logconcave'} with form proportional to exp(-f(x)) where f(x) is L-smooth
#'       and m-strongly-convex, (d) \code{'exponential'} for the exponential distribution.
#'       The default target distribution is the uniform distribution.}
#'     \item{\code{variance}}{The variance of the multidimensional spherical Gaussian or the
#'       exponential distribution. The default value is 1.}
#'     \item{\code{mode}}{A \eqn{d}-dimensional numerical vector that declares the mode of the
#'       Gaussian distribution. The default choice is the center computed by \code{inner_ball()}.}
#'     \item{\code{bias}}{The bias vector for the exponential distribution. The default vector
#'       is \eqn{c_1 = 1} and \eqn{c_i = 0} for \eqn{i \neq 1}.}
#'     \item{\code{L_}}{Smoothness constant (for logconcave).}
#'     \item{\code{m}}{Strong-convexity constant (for logconcave).}
#'     \item{\code{negative_logprob}}{Negative log-probability (for logconcave).}
#'     \item{\code{negative_logprob_gradient}}{Negative log-probability gradient (for logconcave).}
#'   }
#' @param seed Optional. A fixed seed for the number generator.
#' @param ... Additional arguments passed along to methods.
#'
#' @references Robert L. Smith, "Efficient Monte Carlo Procedures for Generating Points
#'   Uniformly Distributed Over Bounded Regions," Operations Research, 1984.
#' @references B.T. Polyak, E.N. Gryazina, "Billiard walk - a new sampling algorithm for control
#'   and optimization," IFAC Proceedings Volumes, 2014.
#' @references Y. Chen, R. Dwivedi, M. J. Wainwright and B. Yu, "Fast MCMC Sampling Algorithms on
#'   Polytopes," Journal of Machine Learning Research, 2018.
#' @references Lee, Yin Tat, Ruoqi Shen, and Kevin Tian, "Logsmooth Gradient Concentration and
#'   Tighter Runtimes for Metropolized Hamiltonian Monte Carlo," arXiv:2002.04121, 2020.
#' @references Shen, Ruoqi, and Yin Tat Lee, "The randomized midpoint method for log-concave
#'   sampling," Advances in Neural Information Processing Systems, 2019.
#' @references Augustin Chevallier, Sylvain Pion, Frederic Cazals, "Hamiltonian Monte Carlo with
#'   boundary reflections, and application to polytope volume calculations," Research Report
#'   preprint hal-01919855, 2018.
#'
#' @return A \eqn{d \times n} matrix that contains, column-wise, the sampled points from the
#'   convex polytope P.
#' @examples
#' # uniform distribution from the 3d unit cube in H-representation using ball walk
#' P = gen_cube(3, 'H')
#' points = sample_points(P, n = 100, random_walk = list("walk" = "BaW", "walk_length" = 5))
#'
#' # gaussian distribution from the 2d unit simplex in H-representation with variance = 2
#' A = matrix(c(-1,0,0,-1,1,1), ncol=2, nrow=3, byrow=TRUE)
#' b = c(0,0,1)
#' P = Hpolytope(A=A,b=b)
#' points = sample_points(P, n = 100, distribution = list("density" = "gaussian", "variance" = 2))
#'
#' # uniform points from the boundary of a 2-dimensional random H-polytope
#' P = gen_rand_hpoly(2,20)
#' points = sample_points(P, n = 100, random_walk = list("walk" = "BRDHR"))
#'
#' # For sampling from logconcave densities see the examples directory
#'
#' @export
sample_points <- function(P, n, random_walk = NULL, distribution = NULL, seed = NULL, ...) {
    UseMethod("sample_points")
}

#' @rdname sample_points
#' @export
sample_points.default <- function(P, n, random_walk = NULL, distribution = NULL, seed = NULL, ...) {
    sample_points_internal(P, n, random_walk, distribution, seed)
}

#' @rdname sample_points
#' @export
sample_points.Spectrahedron <- function(P, n, random_walk = NULL, distribution = NULL, seed = NULL, ...) {
    if (!is.numeric(n) || length(n) != 1L || is.na(n) || n < 1L || n != as.integer(n)) {
        stop("'n' must be a single positive integer")
    }
    stop("Spectrahedron sampling is not implemented in Rvolesti yet; this call will currently error.")
}
