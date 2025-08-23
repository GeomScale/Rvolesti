# VolEsti (volume computation and sampling library)

# Copyright (c) 2012-2025 Vissarion Fisikopoulos
# Copyright (c) 2018-2025 Apostolos Chalkis
# Copyright (c) 2025-2025 Iva Janković

# Contributed and/or modified by Iva Janković, as part of Google Summer of Code 2025 program.

# Licensed under GNU LGPL.3, see LICENCE file

# Example script for comparing volesti boundary samplers 

# Import required libraries
library(ggplot2)
library(volesti)
library(hitandrun)

# Generate 100D cube
P = gen_cube(100, 'H')

# Sampling time (Boundary Random Direction Hit and Run)
start_time_brdhr <- Sys.time()
points_brdhr= sample_points(P, n = 10000, random_walk = list("walk" = "BRDHR"))
end_time_brdhr <- Sys.time()

# Calculate Effective Sample size (Boundary Random Direction Hit and Run)
brdhr_ess = ess(points_brdhr)
min_ess_brdhr <- min(brdhr_ess)

# Calculate PSRF (Boundary Random Direction Hit and Run)
brdhr_psrfs = psrf_univariate(points_brdhr)
max_psrf_brdhr = max(brdhr_psrfs)
elapsed_time_brdhr <- end_time_brdhr - start_time_brdhr
time_per_ind_brdhr <- elapsed_time_brdhr / min_ess_brdhr

#
# Sampling time (Boundary Coordinate Direction Hit and Run)
start_time_bcdhr <- Sys.time()
points_bcdhr= sample_points(P, n = 10000, random_walk = list("walk" = "BCDHR"))
end_time_bcdhr <- Sys.time()

# Calculate Effective Sample size (Boundary Coordinate Direction Hit and Run)
bcdhr_ess = ess(points_bcdhr)
min_ess_bcdhr <- min(bcdhr_ess)

# Calculate PSRF (Boundary Coordinate Direction Hit and Run)
bcdhr_psrfs = psrf_univariate(points_bcdhr)
max_psrf_bcdhr = max(bcdhr_psrfs)
elapsed_time_bcdhr <- end_time_bcdhr - start_time_bcdhr
time_per_ind_bcdhr <- elapsed_time_bcdhr / min_ess_bcdhr

#

# Sampling time (Shake and Bake)
start_time_sb <- Sys.time()
points_sb = sample_points(P, n = 10000, random_walk = list("walk" = "SB"))
end_time_sb <- Sys.time()

# Calculate Effective Sample size (Shake and Bake)
sb_ess = ess(points_sb)
min_ess_sb <- min(sb_ess)

# Calculate PSRF (Shake and Bake)
sb_psrfs = psrf_univariate(points_sb)
max_psrf_sb = max(sb_psrfs)
elapsed_time_sb <- end_time_sb - start_time_sb
time_per_ind_sb <- elapsed_time_sb / min_ess_sb

#

# Sampling time (Billiard Shake and Bake)
start_time_bsb <- Sys.time()
points_bsb = sample_points(P, n = 10000, random_walk = list("walk" = "BSB"))
end_time_bsb <- Sys.time()

# Calculate Effective Sample size (Billiard Shake and Bake)
bsb_ess = ess(points_bsb)
min_ess_bsb <- min(bsb_ess)

# Calculate PSRF (Billiard Shake and Bake)
bsb_psrfs = psrf_univariate(points_bsb)
max_psrf_bsb = max(bsb_psrfs)
elapsed_time_bsb <- end_time_bsb - start_time_bsb
time_per_ind_bsb <- elapsed_time_bsb / min_ess_bsb

#Final results 

results <- data.frame(
  Method = c("Boundary RDHR", "Boundary CDHR","Shake and Bake", "Billiard Shake and Bake"),
  Min_ESS = c(min_ess_brdhr, min_ess_bcdhr,min_ess_sb, min_ess_bsb),
  Max_PSRF = c(max_psrf_brdhr, max_psrf_bcdhr, max_psrf_sb, max_psrf_bsb),
  Time_per_ind_sample = c(time_per_ind_brdhr, time_per_ind_bcdhr,time_per_ind_sb, time_per_ind_bsb)
)

print(results)
