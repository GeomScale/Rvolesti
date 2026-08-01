```r
# Example Usage
library(volesti)

# Define negative log-probability (function value only!)
neg_log_prob <- function(x) {
    # Standard normal: f(x) = 0.5 * ||x||^2
    0.5 * sum(x^2)
}

# Create polytope
P <- gen_cube(d = 5)

# Sample with autodiff (no gradient needed)
samples <- sample_points(
    P,
    WalkType = "CDHR",
    N = 1000,
    use_autodiff = TRUE,
    grad_function = neg_log_prob
)

print(paste("Sampled", ncol(samples), "points from dimension", nrow(samples), "polytope"))
```
