#!/usr/bin/python3

# ------------------------------------------------------------------------------
# PART 1: Import modules and variables
# ------------------------------------------------------------------------------

import numpy as np
import os
from pathlib import Path 		# Manage file paths

# Get repo root
import git
repo = git.Repo('.', search_parent_directories=True)

# Import simulations parameters
import sys
sys.path.append(repo.working_tree_dir + "/setup/")
from setup import SAMPLING_PARAMETERS

# ------------------------------------------------------------------------------
# PART 2: Define functions for simulations
# ------------------------------------------------------------------------------

def run_simulations(ROUTINE_PARAMETERS, SAMPLING_PARAMETERS):
	
	# Import sampling parameters
	therm_N = SAMPLING_PARAMETERS["thermalization_samples"]
	N = SAMPLING_PARAMETERS["montecarlo_samples"]
	number_of_beta = SAMPLING_PARAMETERS["number_of_beta"]
	
	# Run Markov-Chain Monte Carlo simulations
	for L, beta_min, beta_max in ROUTINE_PARAMETERS:
		Path(repo.working_tree_dir + f"/simulations/data/L={L}").mkdir(exist_ok=True) # exist_ok prevents errors if the folder already exists

		for beta in np.linspace(beta_min, beta_max, number_of_beta):
		
			julia_script_filepath = repo.working_tree_dir + "/simulations/src/ising2D_cluster.jl"
			output_data_filepath = repo.working_tree_dir + f"/simulations/data/L={L}/beta={beta}.txt"
			shell_command = "julia " + julia_script_filepath + f" {beta} {L} {therm_N} {N} " + output_data_filepath	
			os.system(f"{shell_command}")

# ------------------------------------------------------------------------------
# PART 3: Run simulations
# ------------------------------------------------------------------------------

if __name__ == "__main__":
	
	error_message = "Error: no option specified \n\
Use --compute-ranges as a call option to compute from scratch simulations ranges \n\
Use --read-ranges as a call option not to read from file simulations ranges \n"

	if len(sys.argv) != 2:
		raise ValueError(error_message)
	else:
		compute_ranges = sys.argv[1]

		if compute_ranges == "--compute-ranges":

			'''
			If the user calls for ranges calculation, request routine_parameters
			to compute the variables ROUTINE_PARAMETERS.
			'''

			from routine_parameters import ROUTINE_PARAMETERS

		elif compute_ranges == "--read-ranges":

			'''
			If the user calls for ranges readout, unpack the file
			routine_parameters.txt where the already calculated ranges are stored.
			'''
			
			routine_parameters_filepath = repo.working_tree_dir + "/setup/routine_parameters.txt"
			_, left, right = np.loadtxt(routine_parameters_filepath, delimiter=',', unpack=True)
			ROUTINE_PARAMETERS = []
			
			from setup import SIZES	
			for i in range(len(SIZES)):
				ROUTINE_PARAMETERS.append((SIZES[i], left[i], right[i]))

		else:
			raise ValueError(error_message)
		
		# Comment here for debugging
		run_simulations(ROUTINE_PARAMETERS, SAMPLING_PARAMETERS)
		
		# Uncomment here for debugging
		# print(ROUTINE_PARAMETERS)
