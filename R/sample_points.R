sample_points <- function(P, WalkType = "CDHR", walk_step = 1, N = 1000,
                         use_autodiff = FALSE, grad_function = NULL, ...) { 
    
    # Add validation
    if (!is.logical(use_autodiff)) {
        stop("use_autodiff must be TRUE or FALSE")
    }
    
    if (is.null(grad_function)) {
        stop("grad_function must be provided")
    }
    
    # Add routing logic
    if (use_autodiff) {
        result <- .Call("sample_points_autodiff", P, WalkType, walk_step, N, 
                       grad_function, PACKAGE = "volesti")
    } else {
        result <- .Call("sample_points_standard", P, WalkType, walk_step, N,
                       grad_function, PACKAGE = "volesti")
    }
    
    return(result)
}