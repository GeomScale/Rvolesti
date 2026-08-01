context("spectrahedron volume interface")

make_toy_spectrahedron <- function() {
    A0 <- diag(c(2, 2))
    A1 <- diag(c(1, 0))
    Spectrahedron(matrices = list(A0, A1))
}

test_that("verbosity validation fails for invalid values", {
    toy_S <- make_toy_spectrahedron()
    expect_error(volume(toy_S, verbosity = -1L), "verbosity must be a single non-negative integer")
    expect_error(volume(toy_S, verbosity = 3L), "verbosity must be a single non-negative integer")
    expect_error(volume(toy_S, verbosity = c(1L, 2L)), "verbosity must be a single non-negative integer")
    expect_error(volume(toy_S, verbosity = NA_integer_), "verbosity must be a single non-negative integer")
})

test_that("verbosity controls messaging before the backend stop", {
    toy_S <- make_toy_spectrahedron()

    expect_error(
        expect_message(volume(toy_S, verbosity = 0L), NA),
        "Volume for Spectrahedron"
    )

    expect_error(
        expect_message(volume(toy_S, verbosity = 1L), "Spectrahedron volume estimation is not implemented"),
        "Volume for Spectrahedron"
    )

    expect_error(
        expect_message(volume(toy_S, verbosity = 2L), "Spectrahedron dimension guess"),
        "Volume for Spectrahedron"
    )
})
