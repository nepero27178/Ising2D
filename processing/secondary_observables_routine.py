#!/usr/bin/python3

# ------------------------------------------------------------------------------
# PART 1: Import modules and variables
# ------------------------------------------------------------------------------

import numpy as np
import os
from pathlib import Path 		# Manage file paths
from datetime import datetime

# Get repo root
import git
repo = git.Repo('.', search_parent_directories=True)

# Import simulations parameters
import sys
sys.path.append(repo.working_tree_dir + "/setup/")
from setup import TOPOLOGY, SIZES, SAMPLING_PARAMETERS

# Read which simulations have been performed
routine_parameters_filepath = repo.working_tree_dir + "/setup/routine_parameters_" + TOPOLOGY + ".txt"
_, left, right = np.loadtxt(routine_parameters_filepath, delimiter=',', unpack=True)
ROUTINE_PARAMETERS = []

for i in range(len(SIZES)):
	ROUTINE_PARAMETERS.append((SIZES[i], left[i], right[i]))

# ------------------------------------------------------------------------------
# PART 2: Define functions for computation of secondary observables
# ------------------------------------------------------------------------------

def compute_secondary_observables(TOPOLOGY, ROUTINE_PARAMETERS, SAMPLING_PARAMETERS):

	'''
	The function compute_secondary_observables uses the blocked data file in
	directory \"Ising2D/processing/data/<your-topology>/
								data-blocklength=<your-blocklength>/L=<your-L>/\"
	and	runs the julia scripts performing resampling.
	The output is a file, stored in \"Ising2D/analysis/data=[...]/\", which
	contains the observables as columns and an info about the used blocklength.
	'''

	# Import sampling parameters
	number_of_beta = SAMPLING_PARAMETERS["number_of_beta"]
	block_optimal_length = SAMPLING_PARAMETERS["block_optimal_length"]
	fake_samples = SAMPLING_PARAMETERS["resampling_fake_samples"]
	
	print(f"\nComputing Binder's cumulant and magnetic susceptibility on topology: {TOPOLOGY}, using blocked data with block length: {block_optimal_length}\n")

	for L, beta_min, beta_max in ROUTINE_PARAMETERS:

		# Write the header line to the file with the current date and time
		output_data_filepath = repo.working_tree_dir + f"/analysis/data/{TOPOLOGY}/L={L}.txt"
		Path(repo.working_tree_dir + f"/analysis/data/{TOPOLOGY}").mkdir(exist_ok=True)
		
		with open(output_data_filepath, "w") as output_data_file:
		    output_data_file.write(f"# beta, U, e_U, chi e_chi [calculated {datetime.now()} on topology: {TOPOLOGY}, using blocked data with block length: {block_optimal_length}]\n")

		# Run the julia resampling script
		for beta in np.linspace(beta_min, beta_max, number_of_beta):
		
			julia_script_filepath = repo.working_tree_dir + "/processing/src/resampling.jl"												# Run julia script
			input_data_filepath = repo.working_tree_dir + f"/processing/data/{TOPOLOGY}/blocklength={block_optimal_length}/L={L}/beta={beta}.txt"	# Use blocked data
			# input_data_filepath = repo.working_tree_dir + f"/simulations/data/{TOPOLOGY}/L={L}/beta={beta}.txt"									# Use raw data
			shell_command = "julia " + julia_script_filepath + f" {L} {beta} {fake_samples} " + input_data_filepath	+ " " + output_data_filepath
			os.system(f"{shell_command}")
			
	print("Done! Check the files at " + repo.working_tree_dir + f"/analysis/data/{TOPOLOGY}/\n")
	return None
		    
# ------------------------------------------------------------------------------
# PART 3: Run calculations
# ------------------------------------------------------------------------------

if __name__ == "__main__":
	
	# Comment here for debugging
	compute_secondary_observables(TOPOLOGY, ROUTINE_PARAMETERS, SAMPLING_PARAMETERS)
	
	# Uncomment here for debugging
	# print(ROUTINE_PARAMETERS)
