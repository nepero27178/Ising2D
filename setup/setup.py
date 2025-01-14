#!/usr/bin/python3

import numpy as np

# Lattice sizes to be simulated
MIN_SIZE = 20
MAX_SIZE = 20
STEP = 20
SIZES = range(MIN_SIZE, MAX_SIZE + STEP, STEP)

# Maximum deviation from interval center (beta_pc) of the simulation temperature range
MAX_DEVIATION= 0.025

# Sampling parameters
SAMPLING_PARAMETERS = {
	"thermalization_samples" : 1000,							# How many steps to thermalize?
	"montecarlo_samples" : 10000,								# How many output raw data?
	"number_of_beta" : 50,										# How many temperatures per simulation?
	"block_optimal_length" : 100,								# How long is the optimal block?
	"block_trial_lengths" : np.linspace(10,150,15,dtype=int),	# Try different block lengths
"resampling_fake_samples" : 100									# How many times do we resample?
}

# Expected magnetic susceptibility parameters
THEORETICAL_CRITICAL_PARAMETERS = {
	# Physical parameters
	"beta_c" : np.log(3)/4,
	"nu" : 1,
	"gamma" : 7/4,
	# Numerical parameters
	"x_0" : -0.1,
	"y_0" : 0.1,
	"c_0" : 0
}

# TODO: make this file imported everywhere
