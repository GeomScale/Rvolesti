#' Apply rounding to a convex polytope (H-polytope, V-polytope or a zonotope)
#'
#' Given a convex H or V polytope or a zonotope as input this function brings the polytope in rounded position based on minimum volume enclosing ellipsoid of a pointset.
#'
#' @param P A convex polytope. It is an object from class (a) Hpolytope or (b) Vpolytope or (c) Zonotope.
#' @param settings Optional. A list of settings.
#' \itemize{
#' \item{\code{method} }{ The method to use for rounding, a) \code{'min_ellipsoid'} for the method based on mimimmum volume enclosing ellipsoid of a dataset, b) \code{'max_ellipsoid'} for the method based on maximum volume enclosed ellipsoid, (c) \code{'isotropy'} for the method based on svd decomposition. The default method is \code{'mee'} for all the representations.}
#' \item{\code{seed} }{ Optional. A fixed seed for the number generator.}
#' }
#'
#' @return A list with 4 elements: (a) a polytope of the same class as the input polytope class and (b) the element "T" which is the matrix of the inverse linear transformation that is applied on the input polytope, (c)  the element "shift" which is the opposite vector of that which has shifted the input polytope, (d) the element "round_value" which is the determinant of the square matrix of the linear transformation that is applied on the input polytope.
#'
#' @references \cite{I.Z.Emiris and V. Fisikopoulos,
#' \dQuote{Practical polytope volume approximation,} \emph{ACM Trans. Math. Soft.,} 2018.},
#' @references \cite{Michael J. Todd and E. Alper Yildirim,
#' \dQuote{On Khachiyan’s Algorithm for the Computation of Minimum Volume Enclosing Ellipsoids,} \emph{Discrete Applied Mathematics,} 2007.}
#' @references \cite{B. Cousins and S. Vempala,
#' \dQuote{A practical volume algorithm,} \emph{Math. Prog. Comp.,} 2016.},
#' @references \cite{Yin Zhang and Liyan Gao,
#' \dQuote{On Numerical Solution of the Maximum Volume Ellipsoid Problem,} \emph{SIAM Journal on Optimization,} 2003.},
#'
#'
#' @examples
#' # round a 5d skinny cube
#' P = gen_skinny_cube(5)
#' listHpoly = round_polytope(P)
#'
#' # round a V-polytope (3d unit cube)
#' P = gen_cube(3, 'V')
#' ListVpoly = round_polytope(P)
#'
#' # round a 2-dimensional zonotope defined by 6 generators
#' Z = gen_rand_zonotope(2,6)
#' ListZono = round_polytope(Z)
#' @export
round_polytope <- function(P, settings = list()){

  seed = NULL
  if (!is.null(settings$seed)) {
    seed = settings$seed
  }

  ret_list = rounding(P, settings$method, seed)

  #get the matrix that describes the polytope
  Mat = ret_list$Mat

  # first column is the vector b
  b = Mat[, 1]

  # remove first column
  A = Mat[, -c(1), drop = FALSE]

  type = P@type
  if (type == 'Vpolytope') {
    PP = list("P" = Vpolytope(V = A), "T" = ret_list$T, "shift" = ret_list$shift, "round_value" = ret_list$round_value)
  }else if (type == 'Zonotope') {
    PP = list("P" = Zonotope(G = A), "T" = ret_list$T, "shift" = ret_list$shift, "round_value" = ret_list$round_value)
  } else {
    #if (dim(P$Aeq)[1] > 0){
    #  PP = list("P" = Hpolytope(A = A, b = b), "T" = ret_list$T, "shift" = ret_list$shift, "round_value" = ret_list$round_value, "N" = ret_list$N, "N_shift" = ret_list$N_shift, "svd_prod" = ret_list$svd_prod)
    PP = list("P" = Hpolytope(A = A, b = b), "T" = ret_list$T, "shift" = ret_list$shift, "round_value" = ret_list$round_value)
  }
  return(PP)
}
