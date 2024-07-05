# VolEsti (volume computation and sampling library)

# Copyright (c) 2012-2020 Vissarion Fisikopoulos
# Copyright (c) 2018-2020 Apostolos Chalkis
# Copyright (c) 2020-2020 Marios Papachristou

# Contributed and/or modified by Atrayee Samanta, as part of Google Summer of Code 2024 program.

# Licensed under GNU LGPL.3, see LICENCE file

# Example for sampling correlation matrices

#Import required library
library(volesti)

#Input parameters:
# n: The dimension of the correlation matrix
# num_matrices: The number of matrices to generate (default is 1000)
# walk_length: The length of the random walk (default is 1)
# nburns: The number of burn-in steps for the random walk (default is 0)
# validate: Whether to validate the sampled matrices (default is FALSE)

result <- uniform_sample_correlation_matrices(n = 5, num_matrices = 10, walk_length = 2, nburns = 0, validate = TRUE)

# result now contains the list of sampled correlation matrices