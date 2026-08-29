context("Simplex-ball R interface")


.sb_two_component_tetrahedron <- function() {
  list(
    A = rbind(
      c(1, 10, 10), c(1, -10, -10),
      c(-1, 10, -10), c(-1, -10, 10)
    ),
    b = rep(2, 4),
    x0 = rep(0, 3),
    V = rbind(
      c(-2, -2, 2, 2),
      c(-0.2, 0.2, -0.2, 0.2),
      c(-0.2, 0.2, 0.2, -0.2)
    )
  )
}


.sb_circle_triangle <- function() {
  vertical <- 3 * sqrt(3) / 2
  list(
    A = rbind(
      c(-1, 0),
      c(0.5, sqrt(3) / 2),
      c(0.5, -sqrt(3) / 2)
    ),
    b = rep(1.5, 3),
    x0 = c(0, 0),
    V = cbind(c(3, 0), c(-1.5, -vertical), c(-1.5, vertical))
  )
}


# The derived vertices are (1.2, 0), (-2, -1), and (-2, 1). The unit
# circle intersects the triangle in a short right arc and a long left arc.
.sb_unequal_arc_triangle <- function() {
  height <- 1
  span <- 3.2
  right_vertex <- 1.2
  amplitude <- sqrt(height^2 + span^2)
  q <- asin(height * right_vertex / amplitude)
  phi <- atan(height / span)

  list(
    A = rbind(c(-1, 0), c(height, span), c(height, -span)),
    b = c(2, height * right_vertex, height * right_vertex),
    x0 = c(0, 0),
    V = cbind(c(right_vertex, 0), c(-2, -1), c(-2, 1)),
    weights = c((q - phi) / (2 * q), (q + phi) / (2 * q))
  )
}


.sb_fast_sample <- function(fixture, walk = "GCW", seed = 901L,
                            N = 8L, component_weights = NULL, ...) {
  sample_simplex_ball(
    fixture$A, fixture$b, fixture$x0,
    N = N,
    walk = walk,
    walk_length = 2L,
    burnin = 2L,
    psrf_target = 100,
    ess_target = 1e-6,
    n_chains = 2L,
    check_interval = 4L,
    max_iterations = 8L,
    seed = seed,
    regcw_tau = 1,
    max_reflections = 20L,
    component_weights = component_weights,
    ...
  )
}


.expect_sb_samples <- function(result, fixture, N, tolerance = 1e-7) {
  expect_s3_class(result, "simplex_ball_sample")
  expect_true(is.matrix(result$samples))
  expect_identical(dim(result$samples), c(ncol(fixture$A), as.integer(N)))
  expect_true(all(is.finite(result$samples)))

  centered <- sweep(result$samples, 1, fixture$x0, "-")
  expect_equal(
    sqrt(colSums(centered^2)),
    rep(1, N),
    tolerance = tolerance
  )
  expect_lte(max(fixture$A %*% result$samples - fixture$b), tolerance)
}


test_that("public Task 5 signatures and defaults remain stable", {
  sampler_formals <- formals(sample_simplex_ball)
  expected_sampler_args <- c(
    "A", "b", "x0", "N", "walk", "walk_length", "burnin",
    "psrf_target", "ess_target", "n_chains", "check_interval",
    "max_iterations", "seed", "regcw_tau", "max_reflections",
    "component_weights", "volume_draws", "max_volume_attempts", "tolerance"
  )
  expect_true(all(expected_sampler_args %in% names(sampler_formals)))
  expect_equal(eval(sampler_formals$N), 1000)
  expect_identical(eval(sampler_formals$walk), "GCW")
  expect_equal(eval(sampler_formals$walk_length), 10)
  expect_equal(eval(sampler_formals$burnin), 100)
  expect_equal(eval(sampler_formals$psrf_target), 1.1)
  expect_equal(eval(sampler_formals$ess_target), 100)
  expect_equal(eval(sampler_formals$n_chains), 3)

  expect_identical(
    names(formals(find_components)),
    c("V", "x0", "radius", "tolerance")
  )
  expect_identical(
    names(formals(find_starting_point)),
    c("V", "component", "A", "b", "x0", "interior_point",
      "radius", "tolerance")
  )
  expect_identical(
    names(formals(ellipsoid_to_simplex_ball)),
    c("covariance", "target_volatility", "tolerance")
  )
})


test_that("vertices are derived and known components use one-based labels", {
  fixture <- .sb_two_component_tetrahedron()

  explicit <- find_components(fixture$V, fixture$x0)
  derived_vertices <- volesti:::.simplex_ball_derive_vertices(
    fixture$A, fixture$b
  )
  derived <- find_components(derived_vertices, fixture$x0)

  expect_s3_class(explicit, "simplex_ball_components")
  expect_equal(derived_vertices, fixture$V, tolerance = 1e-12)
  expect_identical(explicit$components, list(c(1L, 2L), c(3L, 4L)))
  expect_identical(derived$components, explicit$components)
})


test_that("one feasible sphere starting point is returned per component", {
  fixture <- .sb_two_component_tetrahedron()
  components <- find_components(fixture$V, fixture$x0)$components

  starts <- lapply(
    components,
    function(component) {
      find_starting_point(
        fixture$V, component, fixture$A, fixture$b, fixture$x0
      )
    }
  )

  expect_length(starts, 2L)
  for (start in starts) {
    expect_type(start, "double")
    expect_length(start, 3L)
    expect_true(all(is.finite(start)))
    expect_equal(sqrt(sum((start - fixture$x0)^2)), 1, tolerance = 1e-8)
    expect_lte(max(fixture$A %*% start - fixture$b), 1e-8)
  }
})


for (simplex_ball_walk in c("GCW", "ReGCW")) {
  local({
    walk <- simplex_ball_walk

    test_that(paste(walk, "returns finite feasible exact-N samples"), {
      fixture <- .sb_circle_triangle()
      result <- suppressWarnings(.sb_fast_sample(
        fixture,
        walk = walk,
        seed = if (walk == "GCW") 1101L else 1102L
      ))

      .expect_sb_samples(result, fixture, 8L)
      expect_equal(sum(result$component_weights), 1, tolerance = 1e-12)
      expect_equal(sum(result$component_sample_counts), 8)
    })
  })
}


test_that("an explicit seed makes sampling reproducible", {
  fixture <- .sb_circle_triangle()
  first <- suppressWarnings(.sb_fast_sample(fixture, seed = 777L))
  second <- suppressWarnings(.sb_fast_sample(fixture, seed = 777L))
  different <- suppressWarnings(.sb_fast_sample(fixture, seed = 778L))

  expect_identical(first$samples, second$samples)
  expect_identical(first$component_sample_counts, second$component_sample_counts)
  expect_identical(first$component_weights, second$component_weights)
  expect_false(isTRUE(all.equal(first$samples, different$samples)))
})


test_that("native diagnostics, history, and cap warnings are propagated", {
  fixture <- .sb_circle_triangle()
  ordinary <- suppressWarnings(.sb_fast_sample(fixture, seed = 1201L))

  expect_true(all(c(
    "samples", "components", "starting_points", "component_weights",
    "component_sample_counts", "component_results", "diagnostics",
    "diagnostic_history", "converged", "warnings", "call"
  ) %in% names(ordinary)))

  component <- ordinary$component_results[[1L]]
  expect_true(all(c(
    "samples", "chains", "marginal_psrf", "multivariate_psrf",
    "effective_sample_size", "per_chain_ess", "diagnostic_history",
    "diagnostic_checks", "draws_per_chain", "converged", "warning"
  ) %in% names(component)))
  expect_length(component$marginal_psrf, 2L)
  expect_length(component$effective_sample_size, 2L)
  expect_gte(component$diagnostic_checks, 1L)
  expect_length(component$diagnostic_history, component$diagnostic_checks)
  expect_lte(component$draws_per_chain, 8L)

  expect_warning(
    capped <- sample_simplex_ball(
      fixture$A, fixture$b, fixture$x0,
      N = 8L,
      walk = "GCW",
      walk_length = 1L,
      burnin = 1L,
      psrf_target = 1.000001,
      ess_target = 1e9,
      n_chains = 2L,
      check_interval = 4L,
      max_iterations = 4L,
      seed = 1301L
    ),
    "maximum|iteration|converg"
  )
  .expect_sb_samples(capped, fixture, 8L)
  expect_false(capped$converged)
  expect_true(length(capped$warnings) >= 1L)
  expect_false(capped$component_results[[1L]]$converged)
  expect_match(
    capped$component_results[[1L]]$warning,
    "maximum|iteration|converg"
  )
})


test_that("provided and automatic component weights are normalized", {
  fixture <- .sb_unequal_arc_triangle()

  provided <- suppressWarnings(.sb_fast_sample(
    fixture, seed = 1401L, component_weights = c(1, 3)
  ))
  .expect_sb_samples(provided, fixture, 8L)
  expect_identical(provided$components, list(1L, c(2L, 3L)))
  expect_equal(
    unname(provided$component_weights),
    c(0.25, 0.75),
    tolerance = 1e-12
  )
  expect_equal(sum(provided$component_sample_counts), 8)
  expect_true(provided$native_sampling$parallel_components)
  expect_identical(
    provided$native_sampling$active_components,
    c(TRUE, TRUE)
  )

  automatic <- suppressWarnings(.sb_fast_sample(
    fixture, seed = 1402L, volume_draws = 1000L
  ))
  .expect_sb_samples(automatic, fixture, 8L)
  expect_true(all(is.finite(automatic$component_weights)))
  expect_true(all(automatic$component_weights > 0))
  expect_equal(sum(automatic$component_weights), 1, tolerance = 1e-12)
  expect_lt(automatic$component_weights[[1L]], automatic$component_weights[[2L]])
  expect_lt(
    max(abs(unname(automatic$component_weights) - fixture$weights)),
    0.12
  )
})


test_that("zero-weight components are skipped without affecting convergence", {
  fixture <- .sb_unequal_arc_triangle()
  result <- suppressWarnings(.sb_fast_sample(
    fixture,
    seed = 1403L,
    component_weights = c(1, 0)
  ))

  .expect_sb_samples(result, fixture, 8L)
  expect_identical(unname(result$component_sample_counts), c(8L, 0L))
  expect_true(all(result$sample_components == 1L))
  expect_length(result$component_results[[2L]]$chains, 0L)
  expect_identical(result$component_results[[2L]]$draws_per_chain, 0L)
  expect_true(is.na(result$component_converged[[2L]]))
  expect_identical(
    result$native_sampling$active_components,
    c(TRUE, FALSE)
  )
  expect_true(result$converged)
})


test_that("sampling capacity is checked per allocated component", {
  fixture <- .sb_unequal_arc_triangle()

  result <- suppressWarnings(sample_simplex_ball(
    fixture$A, fixture$b, fixture$x0,
    N = 16L,
    walk = "GCW",
    walk_length = 1L,
    burnin = 1L,
    psrf_target = 100,
    ess_target = 1e-6,
    n_chains = 2L,
    check_interval = 4L,
    max_iterations = 4L,
    seed = 1404L,
    component_weights = c(0.5, 0.5)
  ))

  .expect_sb_samples(result, fixture, 16L)
  expect_identical(unname(result$component_sample_counts), c(8L, 8L))
})


test_that("positive-weight unsampled components prevent convergence", {
  fixture <- .sb_unequal_arc_triangle()

  expect_warning(
    result <- sample_simplex_ball(
      fixture$A, fixture$b, fixture$x0,
      N = 1L,
      walk = "GCW",
      walk_length = 1L,
      burnin = 1L,
      psrf_target = 100,
      ess_target = 1e-6,
      n_chains = 2L,
      check_interval = 4L,
      max_iterations = 4L,
      seed = 1405L,
      component_weights = c(0.999, 0.001)
    ),
    "too small|unsampled|unavailable"
  )

  .expect_sb_samples(result, fixture, 1L)
  expect_identical(unname(result$component_sample_counts), c(1L, 0L))
  expect_identical(result$native_sampling$active_components, c(TRUE, FALSE))
  expect_false(result$converged)
})


test_that("ellipsoid conversion round-trips fixed-volatility portfolios", {
  covariance <- diag(c(0.04, 0.09, 0.16))
  target_volatility <- 0.2
  weights <- matrix(c(9, 16, 16) / 41, ncol = 1L)

  transformed <- ellipsoid_to_simplex_ball(
    covariance,
    target_volatility = target_volatility
  )

  expect_true(all(c(
    "A", "b", "x0", "V", "to_original", "to_canonical"
  ) %in% names(transformed)))
  expect_identical(dim(transformed$A), c(3L, 2L))
  expect_identical(dim(transformed$V), c(2L, 3L))
  expect_true(all(is.finite(c(
    transformed$A, transformed$b, transformed$x0, transformed$V
  ))))

  canonical <- transformed$to_canonical(weights)
  recovered <- transformed$to_original(canonical)

  expect_identical(dim(canonical), c(2L, 1L))
  expect_equal(sqrt(colSums(canonical^2)), 1, tolerance = 1e-8)
  expect_lte(max(transformed$A %*% canonical - transformed$b), 1e-8)
  expect_equal(recovered, weights, tolerance = 1e-8)
  expect_equal(colSums(recovered), 1, tolerance = 1e-10)
  expect_gte(min(recovered), -1e-10)
  expect_equal(
    colSums(recovered * (covariance %*% recovered)),
    target_volatility^2,
    tolerance = 1e-10
  )
})


test_that("malformed geometry and sampler controls fail early", {
  fixture <- .sb_circle_triangle()

  expect_error(
    find_components(fixture$V[-1L, , drop = FALSE], fixture$x0),
    "dimension|coordinates|vertices|V"
  )
  expect_error(
    find_components(fixture$V, c(0, 0, 0)),
    "dimension|simplex|vertices|V"
  )
  bad_V <- fixture$V
  bad_V[1L, 1L] <- Inf
  expect_error(find_components(bad_V, fixture$x0), "finite|Inf|NA")
  expect_error(
    find_components(fixture$V, fixture$x0, radius = 0),
    "radius|positive"
  )
  expect_error(
    find_starting_point(
      fixture$V, 0L, fixture$A, fixture$b, fixture$x0
    ),
    "component|index|one-based"
  )

  base_args <- list(
    A = fixture$A, b = fixture$b, x0 = fixture$x0, N = 8L,
    walk_length = 1L, burnin = 1L, psrf_target = 2, ess_target = 1,
    n_chains = 2L, check_interval = 4L, max_iterations = 4L, seed = 1L
  )
  call_sampler <- function(...) {
    do.call(sample_simplex_ball, utils::modifyList(base_args, list(...)))
  }

  expect_error(call_sampler(N = 0L), "N|positive")
  expect_error(call_sampler(b = fixture$b[-1L]), "b|length")
  bad_A <- fixture$A
  bad_A[1L, 1L] <- Inf
  expect_error(call_sampler(A = bad_A), "A|finite")
  expect_error(call_sampler(walk = "bad"), "walk|GCW|ReGCW")
  expect_error(call_sampler(walk_length = 0L), "walk.length|positive")
  expect_error(call_sampler(burnin = -1L), "burnin|non-negative")
  expect_error(call_sampler(psrf_target = 1), "PSRF|psrf|above one")
  expect_error(call_sampler(ess_target = 0), "ESS|ess|positive")
  expect_error(call_sampler(n_chains = 1L), "chain|two|at least")
  expect_error(call_sampler(check_interval = 3L), "check|four|4")
  expect_error(call_sampler(max_iterations = 3L), "maximum|max_iterations|four|4")
  expect_error(call_sampler(seed = NA_integer_), "seed|finite|NA")
  expect_error(call_sampler(volume_draws = 0L), "volume_draws|integer")
  expect_error(call_sampler(tolerance = 0), "tolerance|positive")
  expect_error(
    call_sampler(walk = "ReGCW", regcw_tau = 0),
    "tau|positive"
  )
  expect_error(
    call_sampler(component_weights = c(0)),
    "weight|positive|sum"
  )
  expect_error(
    call_sampler(component_weights = c(1, 1)),
    "weight|component|length"
  )
})


test_that("invalid covariance matrices and infeasible volatility are rejected", {
  covariance <- diag(c(0.04, 0.09, 0.16))
  nonsymmetric <- covariance
  nonsymmetric[1L, 2L] <- 0.01
  nonfinite <- covariance
  nonfinite[1L, 1L] <- NA_real_

  expect_error(
    ellipsoid_to_simplex_ball(matrix(1, 2, 3), 0.2),
    "square|dimension|covariance"
  )
  expect_error(
    ellipsoid_to_simplex_ball(nonsymmetric, 0.2),
    "symmetric|covariance"
  )
  expect_error(
    ellipsoid_to_simplex_ball(nonfinite, 0.2),
    "finite|NA|covariance"
  )
  expect_error(
    ellipsoid_to_simplex_ball(diag(c(1, 0, 0)), 0.8),
    "positive definite|singular|covariance"
  )
  expect_error(ellipsoid_to_simplex_ball(covariance, 0), "volatility|positive")
  expect_error(
    ellipsoid_to_simplex_ball(covariance, 0.1),
    "infeasible|minimum|volatility"
  )
  expect_error(
    ellipsoid_to_simplex_ball(covariance, 0.5),
    "infeasible|maximum|volatility"
  )
  expect_error(
    ellipsoid_to_simplex_ball(diag(c(0.04, 0.09)), 0.25),
    "at least three|dimension"
  )

  # The affine minimum of this correlated two-asset block is below the
  # long-only minimum (0.10), so a 0.09 target must be rejected by the
  # constrained check rather than merely by the affine-hyperplane check.
  correlated <- rbind(
    c(0.01, 0.019, 0),
    c(0.019, 0.04, 0),
    c(0, 0, 0.25)
  )
  constrained_minimum <- ellipsoid_to_simplex_ball(correlated, 0.15)
  expect_gt(
    constrained_minimum$long_only_minimum_volatility,
    constrained_minimum$minimum_volatility
  )
  expect_error(
    ellipsoid_to_simplex_ball(correlated, 0.09),
    "long-only minimum|minimum volatility"
  )
  expect_error(
    ellipsoid_to_simplex_ball(correlated, 0.5),
    "long-only maximum|maximum volatility"
  )
})


test_that("simplex-ball print, summary, and plot methods work", {
  fixture <- .sb_circle_triangle()
  result <- suppressWarnings(.sb_fast_sample(fixture, seed = 1501L))

  printed <- capture.output(print(result))
  expect_true(length(printed) > 0L)
  expect_match(paste(printed, collapse = " "), "component|sample|PSRF|ESS")

  summary_result <- summary(result)
  expect_true(is.list(summary_result) || is.data.frame(summary_result))
  expect_s3_class(summary_result, "summary_simplex_ball_sample")
  expect_identical(
    c("n_samples", "volume_ratio", "psrf", "ess", "mean", "sd") %in%
      names(summary_result$components),
    rep(TRUE, 6L)
  )
  expect_equal(summary_result$components$n_samples, ncol(result$samples))
  expect_equal(summary_result$components$volume_ratio, 1)
  expect_length(summary_result$components$mean[[1L]], nrow(result$samples))
  expect_length(summary_result$components$sd[[1L]], nrow(result$samples))
  expect_identical(
    getS3method("print", "simplex_ball_sample", optional = TRUE),
    print.simplex_ball_sample
  )
  expect_identical(
    getS3method("summary", "simplex_ball_sample", optional = TRUE),
    summary.simplex_ball_sample
  )
  expect_identical(
    getS3method("print", "summary_simplex_ball_sample", optional = TRUE),
    print.summary_simplex_ball_sample
  )

  plot_file <- tempfile(fileext = ".pdf")
  grDevices::pdf(plot_file)
  on.exit({
    grDevices::dev.off()
    unlink(plot_file)
  }, add = TRUE)
  expect_silent(plot(result))
  expect_true(file.exists(plot_file))
})
