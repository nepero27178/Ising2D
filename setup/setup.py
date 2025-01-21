#!/usr/bin/python3

import numpy as np

# ------------------------------------------------------------------------------
# PART 1: Setup simulations and processing
# ------------------------------------------------------------------------------

# Lattice topology "triangular" / "square" / "hexagonal"
TOPOLOGY = "triangular"

# Lattice sizes to be simulated
MIN_SIZE = 20
MAX_SIZE = 160
STEP = 20
SIZES = range(MIN_SIZE, MAX_SIZE + STEP, STEP)

# TODO CHANGE

# Sampling parameters
SAMPLING_PARAMETERS = {
	"thermalization_samples" : int(1e5),							# How many steps to thermalize?
	"montecarlo_samples" : int(1e6),								# How many output raw data?
	"number_of_beta" : 50,											# How many temperatures per simulation?
	"block_optimal_length" : 100,									# How long is the optimal block?
	"resampling_fake_samples" : 100									# How many times do we resample?
}

STD_DEV_ANALYSIS_SAMPLING_PARAMETERS = {
	"number_of_plots" : 5,											# Thus, realize 4 simulations plus one centered at beta_c			
	"block_trial_lengths" : np.linspace(10,150,15,dtype=int),		# Try different block lengths
}

# ------------------------------------------------------------------------------
# PART 2: Setup analysis
# ------------------------------------------------------------------------------

# Expected magnetic susceptibility parameters
THEORETICAL_CRITICAL_PARAMETERS = {
	# Critical exponents
	"nu" : 1,
	"gamma" : 7/4,
}

# Critical temperature depends on topology

# TODO Complete conditional setting of all parameters? Overkill?

if TOPOLOGY == "square":
	THEORETICAL_CRITICAL_PARAMETERS["beta_c"] = np.log(1 + np.sqrt(2))/2
	THEORETICAL_CRITICAL_PARAMETERS["x_0"] = -0.4
	THEORETICAL_CRITICAL_PARAMETERS["y_0"] = 0.11
	THEORETICAL_CRITICAL_PARAMETERS["c_0"] = 0

elif TOPOLOGY == "triangular":
	THEORETICAL_CRITICAL_PARAMETERS["beta_c"] = np.log(3)/4
	THEORETICAL_CRITICAL_PARAMETERS["x_0"] = -0.25
	THEORETICAL_CRITICAL_PARAMETERS["y_0"] = 0.1
	THEORETICAL_CRITICAL_PARAMETERS["c_0"] = 0
	
	# Maximum simulation deviation from interval center (beta_pc) of the simulation temperature range
	MAX_SIMULATION_DEVIATION = 0.025

elif TOPOLOGY == "hexagonal":
	THEORETICAL_CRITICAL_PARAMETERS["beta_c"] = np.log(2 + np.sqrt(3))/2
	THEORETICAL_CRITICAL_PARAMETERS["x_0"] = -0.45
	THEORETICAL_CRITICAL_PARAMETERS["y_0"] = 0.125
	THEORETICAL_CRITICAL_PARAMETERS["c_0"] = 0
	
	# Maximum simulation deviation from interval center (beta_pc) of the simulation temperature range
	MAX_SIMULATION_DEVIATION = 0.035
	
# Maximum fit deviation from interval center (beta_pc) of the simulation temperature range
MAX_FIT_DEVIATION = MAX_SIMULATION_DEVIATION/2
