#' Estimate the volume of a convex body
#'
#' Glue that routes to the classic polytope volume code plus a Spectrahedron-specific
#' volume() S3 method. Keeping this tiny so the old interface still works.
#'
#' @param P A convex body object such as a polytope or spectrahedron.
#' @param settings Volume settings list for the existing C++ pipeline.
#' @param rounding Optional rounding strategy string.
#' @param ... Extra arguments passed along to the dispatched method.
#'
#' @return A list with numeric entries named volume and log_volume.
#' @export
volume <- function(P, settings = NULL, rounding = NULL, ...) UseMethod("volume")

#' @rdname volume
#' @export
volume.default <- function(P, settings = NULL, rounding = NULL, ...) {
    volume_internal(P, settings, rounding)
}

#' @rdname volume
#' @param verbosity Integer flag: 0 = silent, 1 = coarse prints, 2 = debug spam.
#' @export
volume.Spectrahedron <- function(P, settings = NULL, rounding = NULL, verbosity = 0L, ...) {
    verb_level <- as.integer(verbosity)
    if (length(verb_level) != 1L || is.na(verb_level) || verb_level < 0L || verb_level > 2L) {
        stop("verbosity must be a single non-negative integer: 0, 1, or 2")
    }
    unk_dim <- if ("dim" %in% methods::slotNames(P)) P@dim else as.integer(length(P@matrices) - 1L)
    unk_dim <- as.integer(unk_dim)
    if (is.na(unk_dim) || is.null(unk_dim)) {
        unk_dim <- as.integer(length(P@matrices) - 1L)
    }
    if (verb_level >= 1L) {
        message("Spectrahedron volume estimation is not implemented in Rvolesti yet; this call will currently error.")
    }
    if (verb_level >= 2L) {
        message("Spectrahedron dimension guess: ", unk_dim)
    }
    if (verb_level >= 2L) {
        message("Matrix count: ", length(P@matrices))
    }
    volume_spectrahedra(P@matrices, as.integer(unk_dim), verb_level)
}
