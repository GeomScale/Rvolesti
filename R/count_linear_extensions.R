#' Count Linear Extensions of a Partial Order
#'
#' Approximates the number of linear extensions of a partial order (poset)
#' using volume estimation of the corresponding order polytope.
#' Uses the identity: #LE = n! * Vol(OrderPolytope).
#'
#' @param n Integer. Number of elements in the poset (labeled 0 to n-1).
#' @param relations A matrix with 2 columns, where each row (i,j) represents
#'   the order relation a_i <= a_j. For an antichain (no relations), pass
#'   \code{matrix(ncol=2, nrow=0)}.
#' @param algorithm String. Volume approximation algorithm:
#'   \code{"CB"} (Cooling Balls, default),
#'   \code{"CG"} (Cooling Gaussians), or
#'   \code{"SOB"} (Sequence of Balls).
#' @param error Numeric. Upper bound for the approximation error (default 0.1).
#' @param walk_length Integer. Number of steps per random walk (default 1).
#'
#' @return A list with elements:
#' \describe{
#'   \item{linear_extensions}{Approximate count of linear extensions}
#'   \item{volume}{Approximate volume of the order polytope}
#'   \item{dimension}{Number of elements in the poset}
#' }
#'
#' @references
#' Stanley, R.P. (1986). Two Poset Polytopes.
#' \emph{Discrete & Computational Geometry}, 1, 9-23.
#'
#' @examples
#' \dontrun{
#' # Chain: exactly 1 linear extension
#' count_linear_extensions(4, rbind(c(0,1), c(1,2), c(2,3)))
#'
#' # Antichain: exactly n! = 24 linear extensions
#' count_linear_extensions(4, matrix(ncol=2, nrow=0))
#' }
#'
#' @export
count_linear_extensions <- function(n, relations = matrix(ncol=2, nrow=0),
                                     algorithm = "CB", error = 0.1,
                                     walk_length = 1L) {

    if (!is.matrix(relations)) {
        relations <- matrix(as.integer(relations), ncol = 2, byrow = TRUE)
    }
    storage.mode(relations) <- "integer"

    .Call(`_volesti_count_linear_extensions`,
          as.integer(n), relations, algorithm, error, as.integer(walk_length))
}
