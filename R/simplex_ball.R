.simplex_ball_result_field <- function(result, names, default = NULL) {
  for (name in names) {
    if (!is.null(result[[name]])) {
      return(result[[name]])
    }
  }
  default
}

.simplex_ball_components_result <- function(result, number_of_vertices) {
  if (!is.list(result)) {
    stop("The native component finder returned an invalid result.", call. = FALSE)
  }
  components <- result$components
  if (is.null(components) && length(result) > 0L) {
    components <- result[[1L]]
  }
  if (!is.list(components)) {
    stop("The native component finder did not return a component list.",
         call. = FALSE)
  }

  components <- lapply(seq_along(components), function(index) {
    component <- components[[index]]
    if (!is.numeric(component) || anyNA(component) ||
        any(!is.finite(component)) || any(component != floor(component))) {
      stop("Component vertex indices must be finite integers.", call. = FALSE)
    }
    component <- as.integer(component)
    if (length(component) == 0L || any(component < 1L) ||
        any(component > number_of_vertices) || anyDuplicated(component)) {
      stop("The native component finder returned invalid one-based indices.",
           call. = FALSE)
    }
    component
  })

  if (length(components) > 1L) {
    all_indices <- unlist(components, use.names = FALSE)
    if (anyDuplicated(all_indices)) {
      stop("A simplex vertex was assigned to more than one component.",
           call. = FALSE)
    }
  }

  result$components <- components
  if (!is.null(result$adjacency)) {
    adjacency <- as.matrix(result$adjacency)
    if (!all(dim(adjacency) == number_of_vertices)) {
      stop("The native component adjacency matrix has invalid dimensions.",
           call. = FALSE)
    }
    result$adjacency <- adjacency
  }
  class(result) <- unique(c("simplex_ball_components", class(result)))
  result
}

.simplex_ball_sample_matrix <- function(samples, dimension, name = "samples") {
  if (!is.matrix(samples) || !is.numeric(samples) || any(!is.finite(samples))) {
    stop(sprintf("Native `%s` must be a finite numeric matrix.", name),
         call. = FALSE)
  }
  if (nrow(samples) == dimension) {
    normalized <- samples
  } else if (ncol(samples) == dimension) {
    normalized <- t(samples)
  } else {
    stop(sprintf("Native `%s` has incompatible dimensions.", name),
         call. = FALSE)
  }
  if (ncol(normalized) == 0L) {
    stop(sprintf("Native `%s` contains no retained draws.", name),
         call. = FALSE)
  }
  storage.mode(normalized) <- "double"
  normalized
}

.simplex_ball_component_result <- function(result, dimension) {
  if (!is.list(result)) {
    stop("The native sampler returned an invalid result.", call. = FALSE)
  }
  result$samples <- .simplex_ball_sample_matrix(
    .simplex_ball_result_field(result, "samples"), dimension
  )

  chains <- .simplex_ball_result_field(result, "chains", list())
  if (!is.list(chains)) {
    stop("Native `chains` must be a list of matrices.", call. = FALSE)
  }
  if (length(chains) > 0L) {
    chains <- lapply(seq_along(chains), function(index) {
      .simplex_ball_sample_matrix(
        chains[[index]], dimension, sprintf("chains[[%d]]", index)
      )
    })
  }

  marginal_psrf <- as.numeric(.simplex_ball_result_field(
    result, c("marginal_psrf", "psrf"), rep(NA_real_, dimension)
  ))
  if (length(marginal_psrf) != dimension || anyNA(marginal_psrf)) {
    stop("Native marginal PSRF has invalid dimensions or missing values.",
         call. = FALSE)
  }
  multivariate_psrf <- as.numeric(.simplex_ball_result_field(
    result, c("multivariate_psrf", "mpsrf"), NA_real_
  ))
  if (length(multivariate_psrf) != 1L || is.na(multivariate_psrf)) {
    stop("Native multivariate PSRF must be one non-missing number.",
         call. = FALSE)
  }
  effective_sample_size <- as.numeric(.simplex_ball_result_field(
    result,
    c("effective_sample_size", "ess"),
    rep(NA_real_, dimension)
  ))
  if (length(effective_sample_size) != dimension ||
      anyNA(effective_sample_size)) {
    stop("Native ESS has invalid dimensions or missing values.",
         call. = FALSE)
  }

  converged <- .simplex_ball_result_field(result, "converged", FALSE)
  if (length(converged) != 1L || is.na(converged)) {
    stop("Native `converged` must be one non-missing logical value.",
         call. = FALSE)
  }
  native_warning <- .simplex_ball_result_field(result, "warning", character())
  native_warning <- as.character(native_warning)
  native_warning <- native_warning[!is.na(native_warning) & nzchar(native_warning)]

  history <- .simplex_ball_result_field(
    result, c("diagnostic_history", "history"), list()
  )
  checks <- .simplex_ball_result_field(
    result, c("diagnostic_checks", "checks"), length(history)
  )

  result$chains <- chains
  result$marginal_psrf <- marginal_psrf
  result$multivariate_psrf <- multivariate_psrf
  result$effective_sample_size <- effective_sample_size
  result$ess <- effective_sample_size
  result$per_chain_ess <- .simplex_ball_result_field(
    result, "per_chain_ess", list()
  )
  result$diagnostic_history <- history
  result$history <- history
  result$diagnostic_checks <- checks
  result$checks <- checks
  result$draws_per_chain <- .simplex_ball_result_field(
    result, "draws_per_chain", if (length(chains)) ncol(chains[[1L]]) else NA_integer_
  )
  result$converged <- isTRUE(converged)
  result$warning <- native_warning
  result
}

.simplex_ball_empty_component_result <- function(dimension) {
  list(
    samples = matrix(numeric(), nrow = dimension, ncol = 0L),
    chains = list(),
    marginal_psrf = rep(NA_real_, dimension),
    multivariate_psrf = NA_real_,
    effective_sample_size = rep(NA_real_, dimension),
    ess = rep(NA_real_, dimension),
    per_chain_ess = list(),
    diagnostic_history = list(),
    history = list(),
    diagnostic_checks = 0L,
    checks = 0L,
    draws_per_chain = 0L,
    converged = NA,
    warning = character()
  )
}

.simplex_ball_normalize_weights <- function(weights, number_of_components) {
  if (!is.numeric(weights) || length(weights) != number_of_components ||
      anyNA(weights) || any(!is.finite(weights)) || any(weights < 0)) {
    stop(
      sprintf(
        "`component_weights` must contain %d finite, non-negative values.",
        number_of_components
      ),
      call. = FALSE
    )
  }
  total <- sum(weights)
  if (!is.finite(total) || total <= 0) {
    stop("`component_weights` must have a positive finite sum.",
         call. = FALSE)
  }
  as.numeric(weights / total)
}

.simplex_ball_worst_psrf <- function(result) {
  values <- c(result$marginal_psrf, result$multivariate_psrf)
  if (length(values) == 0L || all(is.na(values))) {
    return(NA_real_)
  }
  max(values, na.rm = TRUE)
}

.simplex_ball_minimum_ess <- function(result) {
  values <- result$effective_sample_size
  if (length(values) == 0L || all(is.na(values))) {
    return(NA_real_)
  }
  min(values, na.rm = TRUE)
}

.simplex_ball_history_psrf <- function(result) {
  history <- result$diagnostic_history
  if (is.data.frame(history)) {
    for (name in c("multivariate_psrf", "mpsrf", "psrf")) {
      if (!is.null(history[[name]])) {
        return(as.numeric(history[[name]]))
      }
    }
  }
  if (is.numeric(history) && is.null(dim(history))) {
    return(as.numeric(history))
  }
  if (is.list(history) && length(history) > 0L) {
    values <- vapply(history, function(entry) {
      if (is.list(entry)) {
        multivariate <- .simplex_ball_result_field(
          entry, c("multivariate_psrf", "mpsrf"), NA_real_
        )
        marginal <- .simplex_ball_result_field(
          entry, c("marginal_psrf", "psrf"), numeric()
        )
        candidates <- c(as.numeric(multivariate), as.numeric(marginal))
        if (length(candidates) && !all(is.na(candidates))) {
          return(max(candidates, na.rm = TRUE))
        }
      }
      NA_real_
    }, numeric(1L))
    if (!all(is.na(values))) {
      return(values)
    }
  }
  .simplex_ball_worst_psrf(result)
}

.simplex_ball_format_number <- function(value, digits = 4L) {
  if (length(value) != 1L || is.na(value)) {
    return("NA")
  }
  if (!is.finite(value)) {
    return(as.character(value))
  }
  formatC(value, digits = digits, format = "fg", flag = "#")
}

#' Find connected components of a simplex--sphere intersection
#'
#' Builds the simplex 1-skeleton after removing edges that intersect the ball
#' and vertices strictly inside it. Returned vertex indices follow R's
#' one-based convention.
#'
#' @param V A numeric matrix containing simplex vertices in columns
#'   (`dimension` by `dimension + 1`). Its transpose is also accepted.
#' @param x0 The ball centre.
#' @param radius A strictly positive ball radius. The sampling API uses the
#'   canonical unit radius.
#' @param tolerance A strictly positive numerical tolerance.
#'
#' @return An object of class `simplex_ball_components`, containing
#'   `components`, `adjacency`, `active`, and `active_vertices` when supplied
#'   by the native implementation.
#'
#' @examples
#' V <- rbind(
#'   c(-2, -2,  2,  2),
#'   c(-0.2, 0.2, -0.2, 0.2),
#'   c(-0.2, 0.2,  0.2, -0.2)
#' )
#' components <- find_components(V, x0 = c(0, 0, 0))
#' length(components$components)
#'
#' @export
find_components <- function(V, x0, radius = 1, tolerance = 1e-10) {
  x0 <- .simplex_ball_numeric_vector(x0, "x0")
  if (length(x0) < 2L) {
    stop("`x0` must have dimension at least two.", call. = FALSE)
  }
  radius <- .simplex_ball_positive_scalar(radius, "radius")
  tolerance <- .simplex_ball_positive_scalar(tolerance, "tolerance")
  vertices <- .simplex_ball_normalize_vertices(V, length(x0))
  native <- .simplex_ball_native_function("cpp_simplex_ball_components")
  result <- native(vertices, x0, radius, tolerance)
  .simplex_ball_components_result(result, ncol(vertices))
}

#' Find a feasible starting point on one simplex--sphere component
#'
#' Intersects rays from an interior point toward the component vertices with
#' the sphere and verifies every simplex inequality.
#'
#' @param V A numeric matrix containing simplex vertices in columns. Its
#'   transpose is also accepted.
#' @param component One non-empty vector of one-based simplex vertex indices.
#' @param A A finite constraint matrix.
#' @param b A finite right-hand-side vector defining \eqn{A x \leq b}.
#' @param x0 The ball centre.
#' @param interior_point An optional point strictly inside the simplex--ball
#'   body. If `NULL`, the native implementation computes an approximate
#'   Chebyshev centre.
#' @inheritParams find_components
#'
#' @return A numeric vector on the requested spherical component.
#'
#' @examples
#' V <- rbind(
#'   c(-2, -2,  2,  2),
#'   c(-0.2, 0.2, -0.2, 0.2),
#'   c(-0.2, 0.2,  0.2, -0.2)
#' )
#' A <- rbind(c(1, 10, 10), c(1, -10, -10),
#'            c(-1, 10, -10), c(-1, -10, 10))
#' b <- rep(2, 4)
#' component_data <- find_components(V, c(0, 0, 0))
#' start <- find_starting_point(
#'   V, component_data$components[[1]], A, b, c(0, 0, 0)
#' )
#'
#' @export
find_starting_point <- function(V, component, A, b, x0,
                                interior_point = NULL, radius = 1,
                                tolerance = 1e-10) {
  body <- .simplex_ball_validate_body(A, b, x0)
  vertices <- .simplex_ball_normalize_vertices(V, body$dimension)
  radius <- .simplex_ball_positive_scalar(radius, "radius")
  tolerance <- .simplex_ball_positive_scalar(tolerance, "tolerance")

  if (!is.numeric(component) || length(component) == 0L || anyNA(component) ||
      any(!is.finite(component)) || any(component != floor(component))) {
    stop("`component` must be a non-empty vector of one-based integer indices.",
         call. = FALSE)
  }
  component <- as.integer(component)
  if (any(component < 1L) || any(component > ncol(vertices)) ||
      anyDuplicated(component)) {
    stop("`component` contains invalid or duplicate vertex indices.",
         call. = FALSE)
  }
  if (!is.null(interior_point)) {
    interior_point <- .simplex_ball_numeric_vector(
      interior_point, "interior_point", body$dimension
    )
  }

  native <- .simplex_ball_native_function("cpp_simplex_ball_start")
  point <- as.numeric(native(
    vertices, component, body$A, body$b, body$x0,
    interior_point, radius, tolerance
  ))
  if (length(point) != body$dimension || any(!is.finite(point))) {
    stop("The native starting-point finder returned an invalid point.",
         call. = FALSE)
  }

  verification_tolerance <- max(1e-7, 100 * tolerance)
  if (abs(sqrt(sum((point - body$x0)^2)) - radius) >
      verification_tolerance * max(1, radius)) {
    stop("The native starting point does not lie on the requested sphere.",
         call. = FALSE)
  }
  if (max(as.vector(body$A %*% point) - body$b) > verification_tolerance) {
    stop("The native starting point violates a simplex inequality.",
         call. = FALSE)
  }
  point
}

#' Sample from a simplex--ball intersection with convergence diagnostics
#'
#' Finds every connected component, obtains one feasible start per component,
#' runs independent convergence-controlled chains for all active components in
#' parallel, and combines exactly `N` retained samples using component
#' surface-area weights.
#'
#' @param A A finite `m` by `d` constraint matrix. It must be an irredundant
#'   H-representation of a full-dimensional simplex, with `m = d + 1`.
#' @param b A finite vector of length `nrow(A)` defining
#'   \eqn{A x \leq b}.
#' @param x0 A finite vector of length `ncol(A)` giving the unit-ball centre.
#' @param N The exact total number of samples to return.
#' @param walk Either `"GCW"` or `"ReGCW"`.
#' @param walk_length The number of random-walk transitions between retained
#'   draws.
#' @param burnin The number of burn-in transitions per chain.
#' @param psrf_target A finite PSRF threshold strictly greater than one.
#' @param ess_target A finite, strictly positive ESS threshold.
#' @param n_chains The number of independent diagnostic chains per component;
#'   it must be at least two.
#' @param check_interval Number of retained draws per chain between diagnostic
#'   checks. If `NULL`, a value between 4 and 50 is chosen from `N`.
#' @param max_iterations Maximum retained draws per chain. If `NULL`, it is at
#'   least `check_interval`, `ceiling(N / n_chains)`, and twice the requested
#'   ESS divided across the chains.
#' @param seed Optional non-negative integer seed. Independent master streams
#'   are derived for component-weight estimation and the native sampler. The
#'   native stream then deterministically seeds each component and performs
#'   final sample selection.
#' @param regcw_tau Positive trajectory scale used by `"ReGCW"`.
#' @param max_reflections Non-negative reflection cap; zero selects the native
#'   default.
#' @param component_weights Optional non-negative component surface-area
#'   weights. They are normalized internally. When omitted for a
#'   multi-component body, the native component-weight estimator is used.
#' @param volume_draws Number of accepted full-sphere directions requested by
#'   the component-weight estimator.
#' @param max_volume_attempts Maximum sphere directions attempted by the
#'   component-weight estimator. If `NULL`, it defaults to 100 times
#'   `volume_draws` (and never below `volume_draws`).
#' @param tolerance Positive numerical tolerance used to recover vertices,
#'   identify components, verify starting points, and classify volume draws.
#'
#' @return An object of class `simplex_ball_sample`. Its `samples` matrix has
#'   exactly `N` columns. It also contains one-based `components`, starting
#'   points, normalized component weights, exact component sample counts, raw
#'   per-component chains and diagnostics, diagnostic histories, volume
#'   estimator metadata, a convergence flag, and preserved warning messages.
#'
#' @details Diagnostics are computed separately inside each connected
#'   component. Gelman--Rubin diagnostics must not compare chains targeting
#'   different disconnected conditional distributions. Native worker threads
#'   therefore parallelize components, while each component retains its own
#'   Task 4 diagnostic state. If a native sampler
#'   reaches `max_iterations`, its partial draws and diagnostics are retained,
#'   the final exact-size sample is formed from that pool, and an R warning is
#'   issued. Component counts use the deterministic largest-remainder method.
#'
#' @examples
#' \dontrun{
#' result <- sample_simplex_ball(
#'   A, b, x0, N = 1000, walk = "GCW", walk_length = 10,
#'   burnin = 100, psrf_target = 1.1, ess_target = 100,
#'   n_chains = 3, seed = 123
#' )
#' print(result)
#' summary(result)
#' plot(result)
#' }
#'
#' @export
sample_simplex_ball <- function(
    A, b, x0, N = 1000,
    walk = "GCW",
    walk_length = 10, burnin = 100,
    psrf_target = 1.1, ess_target = 100,
    n_chains = 3, check_interval = NULL, max_iterations = NULL,
    seed = NULL, regcw_tau = 1, max_reflections = 0,
    component_weights = NULL, volume_draws = 10000,
    max_volume_attempts = NULL, tolerance = 1e-10) {
  call <- match.call()
  body <- .simplex_ball_validate_body(A, b, x0)
  A <- body$A
  b <- body$b
  x0 <- body$x0
  dimension <- body$dimension

  N <- .simplex_ball_count(N, "N", 1L)
  walk <- match.arg(walk, c("GCW", "ReGCW"))
  walk_length <- .simplex_ball_count(walk_length, "walk_length", 1L)
  burnin <- .simplex_ball_count(burnin, "burnin", 0L)
  n_chains <- .simplex_ball_count(n_chains, "n_chains", 2L)
  max_reflections <- .simplex_ball_count(
    max_reflections, "max_reflections", 0L
  )
  volume_draws <- .simplex_ball_count(volume_draws, "volume_draws", 1L)
  seed <- .simplex_ball_seed(seed)
  tolerance <- .simplex_ball_positive_scalar(tolerance, "tolerance")

  psrf_target <- .simplex_ball_positive_scalar(psrf_target, "psrf_target")
  if (psrf_target <= 1) {
    stop("`psrf_target` must be strictly greater than one.", call. = FALSE)
  }
  ess_target <- .simplex_ball_positive_scalar(ess_target, "ess_target")
  if (identical(walk, "ReGCW")) {
    regcw_tau <- .simplex_ball_positive_scalar(regcw_tau, "regcw_tau")
  } else {
    # The GCW kernel does not use the reflective trajectory scale.
    regcw_tau <- 1
  }

  if (is.null(check_interval)) {
    check_interval <- max(4L, min(50L, as.integer(ceiling(N / n_chains))))
  } else {
    check_interval <- .simplex_ball_count(
      check_interval, "check_interval", 4L
    )
  }
  if (is.null(max_iterations)) {
    # Total ESS cannot exceed the total number of retained draws.  Use twice
    # the requested ESS as a practical default allowance for autocorrelation,
    # while still retaining enough points to construct the exact-N result.
    default_iterations <- max(
      4,
      as.double(check_interval),
      ceiling(as.double(N) / n_chains),
      ceiling(2 * ess_target / n_chains)
    )
    if (!is.finite(default_iterations) ||
        default_iterations > .Machine$integer.max) {
      stop(
        "`ess_target` is too large for the native iteration counter.",
        call. = FALSE
      )
    }
    max_iterations <- as.integer(default_iterations)
  } else {
    max_iterations <- .simplex_ball_count(
      max_iterations, "max_iterations", 4L
    )
  }
  if (is.null(max_volume_attempts)) {
    attempts <- min(
      .Machine$integer.max,
      max(as.double(volume_draws), 100 * as.double(volume_draws))
    )
    max_volume_attempts <- as.integer(attempts)
  } else {
    max_volume_attempts <- .simplex_ball_count(
      max_volume_attempts, "max_volume_attempts", volume_draws
    )
  }
  if (max_volume_attempts < volume_draws) {
    stop("`max_volume_attempts` must not be smaller than `volume_draws`.",
         call. = FALSE)
  }

  vertices <- .simplex_ball_derive_vertices(A, b, tolerance)
  component_data <- find_components(vertices, x0, 1, tolerance)
  components <- component_data$components
  number_of_components <- length(components)
  if (number_of_components == 0L) {
    stop("The simplex does not intersect the unit sphere.", call. = FALSE)
  }

  starting_points <- lapply(components, function(component) {
    find_starting_point(
      vertices, component, A, b, x0,
      interior_point = NULL, radius = 1, tolerance = tolerance
    )
  })
  starting_point_matrix <- do.call(cbind, starting_points)
  if (!is.matrix(starting_point_matrix)) {
    starting_point_matrix <- matrix(starting_point_matrix, ncol = 1L)
  }

  seeds <- .simplex_ball_seed_stream(seed, 2L)
  volume_seed <- seeds[[1L]]
  sampler_seed <- seeds[[2L]]
  accumulated_warnings <- character()

  if (!is.null(component_weights)) {
    normalized_weights <- .simplex_ball_normalize_weights(
      component_weights, number_of_components
    )
    volume_estimation <- list(
      method = "provided",
      weights = normalized_weights,
      warning = character()
    )
  } else if (number_of_components == 1L) {
    normalized_weights <- 1
    volume_estimation <- list(
      method = "single_component",
      weights = normalized_weights,
      warning = character()
    )
  } else {
    weight_native <- .simplex_ball_native_function(
      "cpp_simplex_ball_component_weights", optional = TRUE
    )
    if (is.null(weight_native)) {
      stop(
        paste0(
          "Multiple components were found, but the native component-weight ",
          "estimator is unavailable. Supply `component_weights`."
        ),
        call. = FALSE
      )
    }
    volume_estimation <- weight_native(
      A, b, vertices, components, starting_point_matrix, x0,
      volume_draws, max_volume_attempts, volume_seed, tolerance
    )
    if (!is.list(volume_estimation) || is.null(volume_estimation$weights)) {
      stop("The native component-weight estimator returned an invalid result.",
           call. = FALSE)
    }
    normalized_weights <- .simplex_ball_normalize_weights(
      volume_estimation$weights, number_of_components
    )
    volume_estimation$method <- "native_full_sphere"
    volume_estimation$weights <- normalized_weights
    volume_warning <- as.character(.simplex_ball_result_field(
      volume_estimation, c("warning", "warnings"), character()
    ))
    volume_warning <- volume_warning[
      !is.na(volume_warning) & nzchar(volume_warning)
    ]
    if (length(volume_warning)) {
      accumulated_warnings <- c(
        accumulated_warnings,
        paste0("Component-weight estimation: ", volume_warning)
      )
    }
    fallback_assignments <- as.numeric(.simplex_ball_result_field(
      volume_estimation, "fallback_assignments", 0
    ))
    if (length(fallback_assignments) != 1L ||
        !is.finite(fallback_assignments) || fallback_assignments < 0) {
      stop(
        "The native component-weight estimator returned invalid metadata.",
        call. = FALSE
      )
    }
    if (fallback_assignments > 0) {
      accumulated_warnings <- c(
        accumulated_warnings,
        sprintf(
          paste0(
            "Component-weight estimation used nearest-start fallback for ",
            "%d ambiguous sphere draws; increase `volume_draws` and review ",
            "the reported weight standard errors."
          ),
          as.integer(fallback_assignments)
        )
      )
    }
    if (any(normalized_weights == 0)) {
      accumulated_warnings <- c(
        accumulated_warnings,
        paste0(
          "At least one detected component received zero Monte Carlo ",
          "weight; increase `volume_draws` or supply `component_weights`."
        )
      )
    }
  }

  component_names <- paste0("component_", seq_len(number_of_components))
  names(normalized_weights) <- component_names

  sampler_native <- .simplex_ball_native_function("cpp_sample_simplex_ball")
  native_sampling_result <- sampler_native(
    A = A,
    b = b,
    V = vertices,
    starting_points = starting_point_matrix,
    component_weights = normalized_weights,
    center = x0,
    N = N,
    n_chains = n_chains,
    check_interval = check_interval,
    max_iterations = max_iterations,
    walk_length = walk_length,
    burnin = burnin,
    psrf_target = psrf_target,
    ess_target = ess_target,
    walk = walk,
    regcw_tau = regcw_tau,
    max_reflections = max_reflections,
    seed = sampler_seed
  )
  native_results <- .simplex_ball_result_field(
    native_sampling_result, "component_results"
  )
  if (!is.list(native_results) ||
      length(native_results) != number_of_components ||
      !isTRUE(.simplex_ball_result_field(
        native_sampling_result, "parallel_components", FALSE
      ))) {
    stop("The native parallel sampler returned an invalid result.",
         call. = FALSE)
  }

  component_counts <- as.integer(.simplex_ball_result_field(
    native_sampling_result, c("component_counts", "sample_counts")
  ))
  native_weights <- as.numeric(.simplex_ball_result_field(
    native_sampling_result, c("component_weights", "weights")
  ))
  active_components <- as.logical(.simplex_ball_result_field(
    native_sampling_result, "active_components"
  ))
  if (length(component_counts) != number_of_components ||
      anyNA(component_counts) || any(component_counts < 0L) ||
      sum(component_counts) != N ||
      length(native_weights) != number_of_components ||
      anyNA(native_weights) || any(!is.finite(native_weights)) ||
      any(native_weights < 0) ||
      abs(sum(native_weights) - 1) > 1e-10 ||
      length(active_components) != number_of_components ||
      anyNA(active_components) ||
      !identical(active_components, component_counts > 0L)) {
    stop("The native parallel sampler returned invalid allocation metadata.",
         call. = FALSE)
  }
  if (max(abs(native_weights - unname(normalized_weights))) > 1e-10) {
    stop("The native parallel sampler changed the normalized weights.",
         call. = FALSE)
  }
  normalized_weights <- native_weights
  names(normalized_weights) <- names(component_counts) <- component_names

  unsampled_positive_components <-
    normalized_weights > 0 & component_counts == 0L
  if (any(unsampled_positive_components)) {
    accumulated_warnings <- c(
      accumulated_warnings,
      paste0(
        "The requested `N` is too small to represent every positive-weight ",
        "component; convergence diagnostics are unavailable for the ",
        "unsampled components."
      )
    )
  }

  component_results <- lapply(seq_len(number_of_components), function(component) {
    if (component_counts[[component]] == 0L) {
      return(.simplex_ball_empty_component_result(dimension))
    }
    .simplex_ball_component_result(native_results[[component]], dimension)
  })
  names(component_results) <- names(normalized_weights)

  for (component in seq_len(number_of_components)) {
    if (component_counts[[component]] == 0L) {
      next
    }

    component_warning <- component_results[[component]]$warning
    if (!component_results[[component]]$converged &&
        length(component_warning) == 0L) {
      component_warning <- paste0(
        "maximum iterations reached before convergence targets"
      )
    }
    if (length(component_warning)) {
      accumulated_warnings <- c(
        accumulated_warnings,
        paste0("Component ", component, ": ", component_warning)
      )
    }
  }

  samples <- .simplex_ball_sample_matrix(
    .simplex_ball_result_field(native_sampling_result, "samples"),
    dimension
  )
  sample_components <- as.integer(.simplex_ball_result_field(
    native_sampling_result,
    c("sample_components", "sample_component")
  ))
  if (ncol(samples) != N || length(sample_components) != N ||
      anyNA(sample_components) || any(sample_components < 1L) ||
      any(sample_components > number_of_components) ||
      !identical(
        as.integer(tabulate(sample_components, number_of_components)),
        unname(component_counts)
      )) {
    stop("The native parallel sampler returned invalid combined samples.",
         call. = FALSE)
  }

  component_psrf <- vapply(
    component_results, .simplex_ball_worst_psrf, numeric(1L)
  )
  component_ess <- vapply(
    component_results, .simplex_ball_minimum_ess, numeric(1L)
  )
  overall_psrf <- if (all(is.na(component_psrf))) {
    NA_real_
  } else {
    max(component_psrf, na.rm = TRUE)
  }
  total_ess <- if (all(is.na(component_ess))) {
    NA_real_
  } else {
    sum(component_ess, na.rm = TRUE)
  }
  sampled_components <- component_counts > 0L
  converged <- !any(unsampled_positive_components) && all(vapply(
    component_results[sampled_components],
    function(result) isTRUE(result$converged),
    logical(1L)
  ))

  accumulated_warnings <- unique(accumulated_warnings[nzchar(
    accumulated_warnings
  )])
  result <- structure(
    list(
      call = call,
      samples = samples,
      sample_components = sample_components,
      sample_component = sample_components,
      A = A,
      b = b,
      x0 = x0,
      V = vertices,
      vertices = vertices,
      components = components,
      component_data = component_data,
      starting_points = starting_points,
      component_weights = normalized_weights,
      component_sample_counts = component_counts,
      component_counts = component_counts,
      component_results = component_results,
      component_samples = lapply(seq_len(number_of_components), function(i) {
        samples[, sample_components == i, drop = FALSE]
      }),
      chains = lapply(component_results, `[[`, "chains"),
      marginal_psrf = lapply(component_results, `[[`, "marginal_psrf"),
      multivariate_psrf = vapply(
        component_results, `[[`, numeric(1L), "multivariate_psrf"
      ),
      effective_sample_size = lapply(
        component_results, `[[`, "effective_sample_size"
      ),
      diagnostic_history = lapply(
        component_results, `[[`, "diagnostic_history"
      ),
      volume_estimation = volume_estimation,
      native_sampling = list(
        parallel_components = TRUE,
        active_components = active_components
      ),
      psrf = overall_psrf,
      total_ess = total_ess,
      ess = total_ess,
      converged = converged,
      warning = if (length(accumulated_warnings)) {
        paste(accumulated_warnings, collapse = "\n")
      } else {
        ""
      },
      warnings = accumulated_warnings,
      parameters = list(
        N = N,
        walk = walk,
        walk_length = walk_length,
        burnin = burnin,
        psrf_target = psrf_target,
        ess_target = ess_target,
        n_chains = n_chains,
        check_interval = check_interval,
        max_iterations = max_iterations,
        seed = seed,
        regcw_tau = regcw_tau,
        max_reflections = max_reflections,
        volume_draws = volume_draws,
        max_volume_attempts = max_volume_attempts,
        tolerance = tolerance
      )
    ),
    class = "simplex_ball_sample"
  )

  # Stable, descriptive aliases keep the high-level object convenient while
  # retaining the names used by the native Task 4 and legacy R prototypes.
  result$settings <- result$parameters
  result$starting_point_matrix <- starting_point_matrix
  result$weight_estimation <- result$volume_estimation
  result$diagnostics <- result$component_results
  result$component_psrf <- component_psrf
  result$component_ess <- component_ess
  result$component_converged <- vapply(
    seq_along(result$component_results),
    function(component) {
      if (result$component_sample_counts[[component]] == 0L) {
        return(NA)
      }
      isTRUE(result$component_results[[component]]$converged)
    },
    logical(1L)
  )

  if (length(accumulated_warnings)) {
    warning(paste(accumulated_warnings, collapse = "\n"), call. = FALSE)
  }
  result
}

#' Print a simplex--ball sampling result
#'
#' @param x A `simplex_ball_sample` object.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @method print simplex_ball_sample
#' @export
print.simplex_ball_sample <- function(x, ...) {
  number <- length(x$components)
  label <- if (number == 1L) "component" else "components"
  cat(sprintf(
    "%d %s found. PSRF = %s. Total ESS = %s.\n",
    number,
    label,
    .simplex_ball_format_number(x$psrf),
    .simplex_ball_format_number(x$total_ess, digits = 6L)
  ))
  cat(sprintf("Returned samples: %d. Converged: %s.\n",
              ncol(x$samples), if (isTRUE(x$converged)) "yes" else "no"))
  if (length(x$warnings)) {
    cat("Warnings:\n", paste0("- ", x$warnings, collapse = "\n"), "\n",
        sep = "")
  }
  invisible(x)
}

#' Summarize simplex--ball samples by connected component
#'
#' @param object A `simplex_ball_sample` object.
#' @param ... Unused.
#' @return An object of class `summary_simplex_ball_sample`. Its component
#'   table contains sample counts, volume ratios, PSRF, ESS, convergence,
#'   coordinate means, and coordinate standard deviations.
#' @method summary simplex_ball_sample
#' @export
summary.simplex_ball_sample <- function(object, ...) {
  number_of_components <- length(object$components)
  dimension <- nrow(object$samples)
  coordinate_names <- rownames(object$samples)
  if (is.null(coordinate_names)) {
    coordinate_names <- paste0("x", seq_len(dimension))
  }

  means <- vector("list", number_of_components)
  standard_deviations <- vector("list", number_of_components)
  for (component in seq_len(number_of_components)) {
    component_samples <- object$samples[
      , object$sample_components == component, drop = FALSE
    ]
    if (ncol(component_samples) == 0L) {
      component_mean <- component_sd <- rep(NA_real_, dimension)
    } else {
      component_mean <- rowMeans(component_samples)
      component_sd <- if (ncol(component_samples) > 1L) {
        apply(component_samples, 1L, stats::sd)
      } else {
        rep(NA_real_, dimension)
      }
    }
    names(component_mean) <- names(component_sd) <- coordinate_names
    means[[component]] <- component_mean
    standard_deviations[[component]] <- component_sd
  }

  table <- data.frame(
    component = seq_len(number_of_components),
    n_samples = as.integer(object$component_sample_counts),
    volume_ratio = as.numeric(object$component_weights),
    psrf = vapply(object$component_results,
                  .simplex_ball_worst_psrf, numeric(1L)),
    ess = vapply(object$component_results,
                 .simplex_ball_minimum_ess, numeric(1L)),
    converged = vapply(
      seq_len(number_of_components),
      function(component) {
        if (object$component_sample_counts[[component]] == 0L) {
          return(NA)
        }
        isTRUE(object$component_results[[component]]$converged)
      },
      logical(1L)
    ),
    stringsAsFactors = FALSE
  )
  table$mean <- I(means)
  table$sd <- I(standard_deviations)

  structure(
    list(
      call = object$call,
      n_samples = ncol(object$samples),
      n_components = number_of_components,
      psrf = object$psrf,
      total_ess = object$total_ess,
      converged = object$converged,
      components = table,
      warnings = object$warnings
    ),
    class = "summary_simplex_ball_sample"
  )
}

#' Print a simplex--ball sampling summary
#'
#' @param x A `summary_simplex_ball_sample` object.
#' @param ... Additional arguments passed to `print.data.frame()` for the
#'   compact component table.
#' @return `x`, invisibly.
#' @method print summary_simplex_ball_sample
#' @export
print.summary_simplex_ball_sample <- function(x, ...) {
  cat("Simplex--ball sample summary\n")
  cat(sprintf("Samples: %d; components: %d; PSRF: %s; total ESS: %s\n",
              x$n_samples, x$n_components,
              .simplex_ball_format_number(x$psrf),
              .simplex_ball_format_number(x$total_ess, digits = 6L)))
  compact <- x$components[
    , c("component", "n_samples", "volume_ratio", "psrf", "ess",
        "converged"),
    drop = FALSE
  ]
  print.data.frame(compact, row.names = FALSE, ...)
  for (component in seq_len(nrow(x$components))) {
    cat(sprintf("Component %d mean: %s\n", component, paste(
      paste0(names(x$components$mean[[component]]), "=",
             format(x$components$mean[[component]], digits = 5L)),
      collapse = ", "
    )))
    cat(sprintf("Component %d sd:   %s\n", component, paste(
      paste0(names(x$components$sd[[component]]), "=",
             format(x$components$sd[[component]], digits = 5L)),
      collapse = ", "
    )))
  }
  if (length(x$warnings)) {
    cat("Warnings:\n", paste0("- ", x$warnings, collapse = "\n"), "\n",
        sep = "")
  }
  invisible(x)
}

#' Plot simplex--ball convergence and component weights
#'
#' Draws the worst PSRF trace for every component and a pie chart of normalized
#' component surface-area weights using base graphics.
#'
#' @param x A `simplex_ball_sample` object.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
plot.simplex_ball_sample <- function(x, ...) {
  number_of_components <- length(x$component_results)
  traces <- lapply(x$component_results, .simplex_ball_history_psrf)
  maximum_checks <- max(1L, max(vapply(traces, length, integer(1L))))
  finite_values <- unlist(traces, use.names = FALSE)
  finite_values <- finite_values[is.finite(finite_values)]
  upper <- max(c(1.1, x$parameters$psrf_target, finite_values), na.rm = TRUE)
  lower <- min(c(1, finite_values), na.rm = TRUE)
  if (!is.finite(lower) || !is.finite(upper) || lower == upper) {
    lower <- 0.95
    upper <- 1.15
  }
  padding <- max(0.02, 0.05 * (upper - lower))
  colours <- grDevices::rainbow(number_of_components)

  old_parameters <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_parameters), add = TRUE)
  graphics::par(mfrow = c(1, 2))
  graphics::plot(
    NA_real_, NA_real_,
    xlim = c(1, maximum_checks),
    ylim = c(lower - padding, upper + padding),
    xlab = "Diagnostic check",
    ylab = "Worst PSRF",
    main = "PSRF trace"
  )
  for (component in seq_len(number_of_components)) {
    trace <- traces[[component]]
    trace[!is.finite(trace)] <- NA_real_
    graphics::lines(seq_along(trace), trace, type = "o",
                    col = colours[[component]], pch = component)
  }
  graphics::abline(h = x$parameters$psrf_target, lty = 2L, col = "firebrick")
  graphics::legend(
    "topright",
    legend = paste("Component", seq_len(number_of_components)),
    col = colours,
    lty = 1L,
    pch = seq_len(number_of_components),
    bty = "n"
  )
  graphics::pie(
    x$component_weights,
    labels = paste0(
      "Component ", seq_len(number_of_components), "\n",
      formatC(100 * x$component_weights, digits = 3L, format = "fg"), "%"
    ),
    col = colours,
    main = "Component volume ratios"
  )
  invisible(x)
}
