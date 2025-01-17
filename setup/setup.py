#!/usr/bin/python3

import numpy as np

# ------------------------------------------------------------------------------
# PART 1: Setup simulations and processing
# ------------------------------------------------------------------------------

# Lattice topology "triangular" / "square" / "hexagonal"
TOPOLOGY = "triangular"

# Lattice sizes to be simulated
MIN_SIZE = 20
MAX_SIZE = 200
STEP = 20
SIZES = range(MIN_SIZE, MAX_SIZE + STEP, STEP)

# Maximum simulation deviation from interval center (beta_pc) of the simulation temperature range
MAX_SIMULATION_DEVIATION = 0.025

# Sampling parameters
SAMPLING_PARAMETERS = {
	"thermalization_samples" : 1000,							# How many steps to thermalize?
	"montecarlo_samples" : 10000,								# How many output raw data?
	"number_of_beta" : 50,										# How many temperatures per simulation?
	"block_optimal_length" : 1,									# How long is the optimal block?
	"block_trial_lengths" : np.linspace(50,150,40,dtype=int),	# Try different block lengths
	"resampling_fake_samples" : 100								# How many times do we resample?
}

# ------------------------------------------------------------------------------
# PART 2: Setup analysis
# ------------------------------------------------------------------------------

# Maximum fit deviation from interval center (beta_pc) of the simulation temperature range
MAX_FIT_DEVIATION = MAX_SIMULATION_DEVIATION/2

# Expected magnetic susceptibility parameters
THEORETICAL_CRITICAL_PARAMETERS = {
	# Critical exponents
	"nu" : 1,
	"gamma" : 7/4,
	# Numerical parameters
	"x_0" : -0.1,
	"y_0" : 0.1,
	"c_0" : 0
}

# Critical temperature depends on topology

# TODO Complete conditional setting of all parameters? Overkill?

if TOPOLOGY == "square":
	THEORETICAL_CRITICAL_PARAMETERS["beta_c"] = np.log(1 + np.sqrt(2))/2

elif TOPOLOGY == "triangular":
	THEORETICAL_CRITICAL_PARAMETERS["beta_c"] = np.log(3)/4

elif TOPOLOGY == "hexagonal":
	THEORETICAL_CRITICAL_PARAMETERS["beta_c"] = np.log(2 + np.sqrt(3))/2
