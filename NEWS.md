# volesti 1.0.0

* This is the first release of volesti R-Package.

# volesti 1.0.1

* Fix some bugs for solaris os.

# volesti 1.0.2

* Remove r-striper to avoid CRAN policy violation.

# volesti 1.0.3

* Fix CRAN warnings.

# volesti 1.1.0

* New volume computation algorithm.
* Billiard walk for uniform sampling.
* Modified exact volume computation function.
* Implementation and evaluation of PCA method for zonotope approximation.
* Boundary sampling.
* Improved functionality for finance applications.
* Improved names for functions and input variables.
* Use exclusively Eigen/BH library for linear algebra.

# volesti 1.1.1

* Fix CRAN warnings about deprecated use of ftime

# volesti 1.1.2

- Improve functionality of R function.
- Use lower case and "_" separated names for the volesti's functions.
- Use R classes for the convex polytopes

# volesti 1.1.2-2

- Remove `loadSdpaFormatFile()` and `readSdpaFormatFile()` functions.

# volesti 1.1.2-3

- Remove unneeded-internal-declaration warning in clang-15

# volesti 1.1.2-4

- Remove uninitialized warning in clang-16 (lp_presolve)

# volesti 1.1.2-6

- Fix UBSAN issues (lp_presolve)

# volesti 1.2.0

- New functions: dinvweibull_with_loc, ess, estimtate_lipschitz_constant, gen_birkhoff, geweke
pinvweibull_with_loc, psrf_multivariate, psrf_univariate, raftery

- New features in sample_points function:
 a) new walks: i) Dikin walk, ii) Vaidya walk, iii) John walk,
  iv) Hamiltonian Monte Carlo (HMC) for general logconcave densities
  v) Underdamped Langevin Dynamics (ULD) for general logconcave densities (using the Randomized Midpoint Method)
  vi) Exact Hamiltonian Monte Carlo with reflections (spherical Gaussian or exponential distribution)
 b) new distributions: i) exponential, ii) general logconcave
