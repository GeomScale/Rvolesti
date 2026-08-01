```r
test_that("autodiff_produces_samples", {
    P <- gen_cube(d=5)
    f <- function(x) sum(x^2)/2
    
    samples <- sample_points(P, use_autodiff=TRUE, grad_function=f, N=50)
    expect_equal(nrow(samples), 5)
    expect_equal(ncol(samples), 50)
})

test_that("autodiff_matches_manual", {
    P <- gen_cube(d=3)
    f <- function(x) sum(x^2)/2
    grad_f <- function(x) -x
    
    set.seed(42)
    s1 <- sample_points(P, use_autodiff=TRUE, grad_function=f, N=50)
    
    set.seed(42)
    s2 <- sample_points(P, use_autodiff=FALSE, grad_function=grad_f, N=50)
    
    expect_true(max(abs(s1 - s2)) < 0.1)
})
```