# Internal validation and conversion helpers used by the simplex--ball API.

.simplex_ball_numeric_matrix <- function(x, name) {
  if (!is.matrix(x) || !is.numeric(x) || length(dim(x)) != 2L) {
    stop(sprintf("`%s` must be a numeric matrix.", name), call. = FALSE)
  }
  if (any(dim(x) == 0L)) {
    stop(sprintf("`%s` must have positive dimensions.", name), call. = FALSE)
  }
  if (any(!is.finite(x))) {
    stop(sprintf("`%s` must contain only finite values.", name), call. = FALSE)
  }
  storage.mode(x) <- "double"
  x
}

.simplex_ball_numeric_vector <- function(x, name, length = NULL) {
  if (!is.numeric(x) || is.matrix(x) || any(!is.finite(x))) {
    stop(sprintf("`%s` must be a finite numeric vector.", name), call. = FALSE)
  }
  if (!is.null(length) && base::length(x) != length) {
    stop(sprintf("`%s` must have length %d.", name, length), call. = FALSE)
  }
  as.numeric(x)
}

.simplex_ball_count <- function(x, name, minimum = 0L) {
  if (length(x) != 1L || !is.numeric(x) || is.na(x) || !is.finite(x) ||
      x != floor(x) || x < minimum || x > .Machine$integer.max) {
    stop(
      sprintf("`%s` must be one integer greater than or equal to %d.",
              name, minimum),
      call. = FALSE
    )
  }
  as.integer(x)
}

.simplex_ball_positive_scalar <- function(x, name, strictly = TRUE) {
  invalid <- length(x) != 1L || !is.numeric(x) || is.na(x) || !is.finite(x)
  if (strictly) {
    invalid <- invalid || x <= 0
  } else {
    invalid <- invalid || x < 0
  }
  if (invalid) {
    qualifier <- if (strictly) "positive" else "non-negative"
    stop(sprintf("`%s` must be one finite %s number.", name, qualifier),
         call. = FALSE)
  }
  as.numeric(x)
}

.simplex_ball_seed <- function(seed) {
  if (is.null(seed)) {
    return(NULL)
  }
  .simplex_ball_count(seed, "seed", 0L)
}

.simplex_ball_validate_body <- function(A, b, x0, minimum_dimension = 2L) {
  A <- .simplex_ball_numeric_matrix(A, "A")
  d <- ncol(A)
  if (d < minimum_dimension) {
    stop(sprintf("The simplex dimension must be at least %d.", minimum_dimension),
         call. = FALSE)
  }
  b <- .simplex_ball_numeric_vector(b, "b", nrow(A))
  x0 <- .simplex_ball_numeric_vector(x0, "x0", d)

  row_norms <- sqrt(rowSums(A * A))
  if (any(!is.finite(row_norms)) || any(row_norms == 0)) {
    stop("Every row of `A` must have a non-zero finite normal.", call. = FALSE)
  }

  list(A = A, b = b, x0 = x0, dimension = d)
}

.simplex_ball_normalize_vertices <- function(V, dimension,
                                             require_simplex = TRUE) {
  V <- .simplex_ball_numeric_matrix(V, "V")

  if (nrow(V) == dimension) {
    vertices <- V
  } else if (ncol(V) == dimension) {
    vertices <- t(V)
  } else {
    stop(
      sprintf("`V` must have %d coordinates in its rows or columns.",
              dimension),
      call. = FALSE
    )
  }

  if (require_simplex && ncol(vertices) != dimension + 1L) {
    stop(
      sprintf("A full-dimensional %d-simplex must have exactly %d vertices.",
              dimension, dimension + 1L),
      call. = FALSE
    )
  }
  vertices
}

.simplex_ball_derive_vertices <- function(A, b, tolerance = 1e-10) {
  d <- ncol(A)
  if (nrow(A) != d + 1L) {
    stop(
      sprintf(paste0(
        "`A` and `b` must be an irredundant H-representation of a ",
        "%d-simplex, with exactly %d inequalities."), d, d + 1L),
      call. = FALSE
    )
  }

  solve_tolerance <- max(tolerance, sqrt(.Machine$double.eps))
  feasibility_tolerance <- max(
    tolerance,
    100 * .Machine$double.eps * max(1, max(abs(b)))
  )
  vertices <- matrix(NA_real_, nrow = d, ncol = d + 1L)

  for (omitted_facet in seq_len(d + 1L)) {
    keep <- setdiff(seq_len(d + 1L), omitted_facet)
    active_A <- A[keep, , drop = FALSE]
    active_b <- b[keep]
    decomposition <- qr(active_A, tol = solve_tolerance)
    if (decomposition$rank != d) {
      stop(
        sprintf("The facets opposite vertex %d are rank deficient.",
                omitted_facet),
        call. = FALSE
      )
    }

    vertex <- tryCatch(
      as.numeric(qr.solve(active_A, active_b, tol = solve_tolerance)),
      error = function(error) NULL
    )
    if (is.null(vertex) || length(vertex) != d || any(!is.finite(vertex))) {
      stop(sprintf("Could not recover simplex vertex %d.", omitted_facet),
           call. = FALSE)
    }

    residual_scale <- max(1, max(abs(active_b)))
    if (max(abs(as.vector(active_A %*% vertex) - active_b)) >
        feasibility_tolerance * residual_scale ||
        max(as.vector(A %*% vertex) - b) > feasibility_tolerance) {
      stop(
        "`A` and `b` do not describe a bounded, irredundant simplex.",
        call. = FALSE
      )
    }
    vertices[, omitted_facet] <- vertex
  }

  distances <- as.matrix(stats::dist(t(vertices)))
  diag(distances) <- Inf
  vertex_scale <- max(1, max(abs(vertices)))
  if (min(distances) <= feasibility_tolerance * vertex_scale) {
    stop("`A` and `b` produce duplicate simplex vertices.", call. = FALSE)
  }

  edge_matrix <- sweep(vertices[, -1L, drop = FALSE], 1L,
                       vertices[, 1L], "-")
  if (qr(edge_matrix, tol = solve_tolerance)$rank != d) {
    stop("The recovered simplex is not full-dimensional.", call. = FALSE)
  }
  vertices
}

.simplex_ball_native_function <- function(name, optional = FALSE) {
  fun <- get0(name, mode = "function", inherits = TRUE)
  if (is.null(fun) && !optional) {
    stop(
      sprintf("The native `%s()` entry point is not available.", name),
      call. = FALSE
    )
  }
  fun
}

.simplex_ball_with_seed <- function(seed, expression) {
  if (is.null(seed)) {
    return(force(expression))
  }

  global <- .GlobalEnv
  had_seed <- exists(".Random.seed", envir = global, inherits = FALSE)
  if (had_seed) {
    previous_seed <- get(".Random.seed", envir = global, inherits = FALSE)
  }
  on.exit({
    if (had_seed) {
      assign(".Random.seed", previous_seed, envir = global)
    } else if (exists(".Random.seed", envir = global, inherits = FALSE)) {
      rm(".Random.seed", envir = global)
    }
  }, add = TRUE)

  set.seed(seed)
  force(expression)
}

.simplex_ball_seed_stream <- function(seed, number) {
  if (number == 0L) {
    return(integer())
  }
  upper <- .Machine$integer.max - 1L
  .simplex_ball_with_seed(
    seed,
    as.integer(sample.int(upper, number, replace = TRUE) - 1L)
  )
}

.simplex_ball_long_only_minimum_variance <- function(
    base_point, equality_basis, reduced_covariance,
    linear_term, base_variance, tolerance) {
  dimension <- ncol(equality_basis)
  objective <- function(reduced_weights) {
    as.numeric(
      base_variance +
        2 * crossprod(linear_term, reduced_weights) +
        crossprod(
          reduced_weights,
          reduced_covariance %*% reduced_weights
        )
    )
  }
  gradient <- function(reduced_weights) {
    as.numeric(2 * (
      linear_term + reduced_covariance %*% reduced_weights
    ))
  }

  # w = base_point + equality_basis %*% y >= 0.  The equal-weight point
  # y = 0 is strictly feasible, which is required by constrOptim's barrier.
  solution <- tryCatch(
    suppressWarnings(stats::constrOptim(
      theta = numeric(dimension),
      f = objective,
      grad = gradient,
      ui = equality_basis,
      ci = -base_point,
      method = "BFGS",
      control = list(
        reltol = max(tolerance, 1e-12),
        maxit = 2000L
      ),
      outer.iterations = 100L,
      outer.eps = max(tolerance, 1e-10)
    )),
    error = function(error) NULL
  )
  if (is.null(solution) ||
      is.null(solution$convergence) ||
      solution$convergence != 0L ||
      !is.finite(solution$value)) {
    stop(
      "Could not compute the long-only minimum portfolio volatility.",
      call. = FALSE
    )
  }

  weights <- as.numeric(
    base_point + equality_basis %*% solution$par
  )
  feasibility_tolerance <- max(100 * tolerance, 1e-8)
  if (min(weights) < -feasibility_tolerance ||
      abs(sum(weights) - 1) > feasibility_tolerance) {
    stop(
      "The long-only minimum-volatility optimization was infeasible.",
      call. = FALSE
    )
  }

  max(0, as.numeric(solution$value))
}

#' Transform a fixed-volatility portfolio problem to simplex--ball form
#'
#' Reduces the fully-invested portfolio hyperplane to dimension \eqn{n-1},
#' completes the square in the covariance quadratic form, and whitens the
#' resulting ellipsoid to the unit sphere. The transformed feasible set is
#' \eqn{\{z: A z \le b, \|z-x0\|=1\}}.
#'
#' @param covariance A finite, symmetric asset covariance matrix. It must be
#'   positive definite on the hyperplane orthogonal to the all-ones vector.
#' @param target_volatility A strictly positive target standard deviation.
#' @param tolerance Positive numerical tolerance for symmetry, rank and
#'   feasibility checks.
#'
#' @return A list with transformed constraints `A`, `b`, ball centre `x0`,
#'   column-wise simplex vertices `V` (also named `vertices`), and
#'   transformation metadata. The `forward`/`to_canonical` function maps
#'   portfolio-weight vectors or column-wise matrices to transformed
#'   coordinates; `inverse`/`to_original` performs the reverse mapping.
#'   `minimum_volatility` is the minimum volatility on the fully-invested
#'   affine hyperplane before imposing long-only constraints;
#'   `long_only_minimum_volatility` and `long_only_maximum_volatility` give
#'   the feasible volatility interval on the portfolio simplex.
#'
#' @details If \eqn{w=w_0+Qy}, where the columns of \eqn{Q} span
#'   \eqn{\{w:1^T w=0\}}, then the reduced quadratic is completed around
#'   \eqn{\mu=-(Q^T\Sigma Q)^{-1}Q^T\Sigma w_0}. The target variance must be
#'   strictly above the resulting affine-hyperplane minimum; otherwise the
#'   unit-sphere representation is empty or degenerate.
#'
#' @examples
#' Sigma <- diag(c(0.04, 0.09, 0.16))
#' transform <- ellipsoid_to_simplex_ball(Sigma, target_volatility = 0.20)
#' weights <- c(9, 16, 16) / 41
#' transformed <- transform$forward(weights)
#' stopifnot(isTRUE(all.equal(transform$inverse(transformed), weights)))
#'
#' @export
ellipsoid_to_simplex_ball <- function(
    covariance,
    target_volatility,
    tolerance = sqrt(.Machine$double.eps)) {
  covariance <- .simplex_ball_numeric_matrix(covariance, "covariance")
  if (nrow(covariance) != ncol(covariance)) {
    stop("`covariance` must be square.", call. = FALSE)
  }
  n_assets <- nrow(covariance)
  if (n_assets < 3L) {
    stop(
      paste0(
        "`covariance` must describe at least three assets so that the ",
        "canonical simplex--ball body has dimension at least two."
      ),
      call. = FALSE
    )
  }
  target_volatility <- .simplex_ball_positive_scalar(
    target_volatility, "target_volatility"
  )
  tolerance <- .simplex_ball_positive_scalar(tolerance, "tolerance")

  covariance_scale <- max(1, max(abs(covariance)))
  symmetry_error <- max(abs(covariance - t(covariance)))
  symmetry_tolerance <- tolerance * covariance_scale
  if (symmetry_error > symmetry_tolerance) {
    stop("`covariance` must be symmetric.", call. = FALSE)
  }
  covariance <- (covariance + t(covariance)) / 2

  covariance_eigenvalues <- eigen(
    covariance,
    symmetric = TRUE,
    only.values = TRUE
  )$values
  if (min(covariance_eigenvalues) < -tolerance * covariance_scale) {
    stop("`covariance` must be positive semidefinite.", call. = FALSE)
  }

  base_point <- rep(1 / n_assets, n_assets)
  full_basis <- qr.Q(qr(matrix(1, nrow = n_assets, ncol = 1L)),
                     complete = TRUE)
  equality_basis <- full_basis[, -1L, drop = FALSE]
  reduced_covariance <- crossprod(
    equality_basis,
    covariance %*% equality_basis
  )
  reduced_covariance <- (reduced_covariance + t(reduced_covariance)) / 2

  decomposition <- eigen(reduced_covariance, symmetric = TRUE)
  eigen_scale <- max(1, max(abs(decomposition$values)))
  eigen_tolerance <- tolerance * eigen_scale
  if (any(decomposition$values <= eigen_tolerance)) {
    stop(
      paste0(
        "`covariance` must be positive definite on the fully-invested ",
        "portfolio hyperplane."
      ),
      call. = FALSE
    )
  }

  linear_term <- as.numeric(crossprod(
    equality_basis,
    covariance %*% base_point
  ))
  reduced_center <- -as.numeric(solve(reduced_covariance, linear_term))
  base_variance <- as.numeric(crossprod(base_point, covariance %*% base_point))
  minimum_variance <- base_variance - as.numeric(crossprod(
    reduced_center,
    reduced_covariance %*% reduced_center
  ))
  minimum_variance <- max(0, minimum_variance)
  target_variance <- target_volatility^2
  radius_squared <- target_variance - minimum_variance
  variance_tolerance <- tolerance *
    max(1, abs(target_variance), abs(minimum_variance))
  if (!is.finite(radius_squared) || radius_squared <= variance_tolerance) {
    stop(
      sprintf(
        paste0(
          "`target_volatility` must be strictly greater than the ",
          "fully-invested minimum volatility (%.10g)."
        ),
        sqrt(minimum_variance)
      ),
      call. = FALSE
    )
  }

  maximum_simplex_variance <- max(diag(covariance))
  long_only_minimum_variance <-
    .simplex_ball_long_only_minimum_variance(
      base_point,
      equality_basis,
      reduced_covariance,
      linear_term,
      base_variance,
      tolerance
    )
  simplex_variance_tolerance <- tolerance * max(
    1,
    abs(target_variance),
    abs(long_only_minimum_variance),
    abs(maximum_simplex_variance)
  )
  if (target_variance <=
      long_only_minimum_variance + simplex_variance_tolerance) {
    stop(
      sprintf(
        paste0(
          "`target_volatility` must be strictly greater than the ",
          "long-only minimum volatility (%.10g)."
        ),
        sqrt(long_only_minimum_variance)
      ),
      call. = FALSE
    )
  }
  if (target_variance >=
      maximum_simplex_variance - simplex_variance_tolerance) {
    stop(
      sprintf(
        paste0(
          "`target_volatility` must be strictly smaller than the ",
          "long-only maximum volatility (%.10g)."
        ),
        sqrt(maximum_simplex_variance)
      ),
      call. = FALSE
    )
  }

  eigenvectors <- decomposition$vectors
  whitener <- eigenvectors %*%
    diag(sqrt(decomposition$values / radius_squared), nrow = n_assets - 1L) %*%
    t(eigenvectors)
  inverse_whitener <- eigenvectors %*%
    diag(sqrt(radius_squared / decomposition$values),
         nrow = n_assets - 1L) %*%
    t(eigenvectors)

  reduced_A <- -equality_basis
  reduced_b <- base_point
  A <- reduced_A %*% inverse_whitener
  b <- reduced_b - as.numeric(reduced_A %*% reduced_center)
  x0 <- rep(0, n_assets - 1L)

  original_vertices <- diag(n_assets)
  centered_vertices <- sweep(original_vertices, 1L, base_point, "-")
  reduced_vertices <- crossprod(equality_basis, centered_vertices)
  V <- whitener %*% sweep(reduced_vertices, 1L, reduced_center, "-")

  forward <- function(weights) {
    vector_input <- is.numeric(weights) && is.null(dim(weights))
    if (vector_input) {
      weights <- matrix(weights, ncol = 1L)
    }
    weights <- .simplex_ball_numeric_matrix(weights, "weights")
    if (nrow(weights) != n_assets) {
      stop(sprintf("`weights` must have %d rows.", n_assets), call. = FALSE)
    }
    reduced <- crossprod(
      equality_basis,
      sweep(weights, 1L, base_point, "-")
    )
    transformed <- whitener %*%
      sweep(reduced, 1L, reduced_center, "-")
    if (vector_input) as.numeric(transformed) else transformed
  }

  inverse <- function(points) {
    vector_input <- is.numeric(points) && is.null(dim(points))
    if (vector_input) {
      points <- matrix(points, ncol = 1L)
    }
    points <- .simplex_ball_numeric_matrix(points, "points")
    if (nrow(points) != n_assets - 1L) {
      stop(sprintf("`points` must have %d rows.", n_assets - 1L),
           call. = FALSE)
    }
    reduced <- sweep(inverse_whitener %*% points, 1L,
                     reduced_center, "+")
    weights <- sweep(equality_basis %*% reduced, 1L, base_point, "+")
    if (vector_input) as.numeric(weights) else weights
  }

  structure(
    list(
      A = A,
      b = as.numeric(b),
      x0 = x0,
      V = V,
      vertices = V,
      forward = forward,
      inverse = inverse,
      to_canonical = forward,
      to_original = inverse,
      covariance = covariance,
      target_volatility = target_volatility,
      target_variance = target_variance,
      minimum_volatility = sqrt(minimum_variance),
      minimum_variance = minimum_variance,
      long_only_minimum_volatility = sqrt(long_only_minimum_variance),
      long_only_minimum_variance = long_only_minimum_variance,
      long_only_maximum_volatility = sqrt(maximum_simplex_variance),
      long_only_maximum_variance = maximum_simplex_variance,
      base_point = base_point,
      equality_basis = equality_basis,
      reduced_center = reduced_center,
      reduced_covariance = reduced_covariance,
      whitener = whitener,
      inverse_whitener = inverse_whitener
    ),
    class = "simplex_ball_transform"
  )
}
