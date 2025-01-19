#!/usr/bin/python3

'''
This script generates very few ready-to-use data to implement std-dev-analysis.
Raw data are not directly used because they refer to different betas, thus
making it impossible to plot data from the same temperature for different
sizes on the same plot.
'''

# ------------------------------------------------------------------------------
# PART 1: Import modules and variables
# ------------------------------------------------------------------------------

import numpy as np
import os
from datetime import datetime	
from pathlib import Path 		# Manage file paths

# Get repo root
import git
repo = git.Repo('.', search_parent_directories=True)

# Import simulations parameters
import sys
sys.path.append(repo.working_tree_dir + "/setup/")
from setup import TOPOLOGY, SIZES, MAX_SIMULATION_DEVIATION, THEORETICAL_CRITICAL_PARAMETERS, SAMPLING_PARAMETERS, STD_DEV_ANALYSIS_SAMPLING_PARAMETERS

# ------------------------------------------------------------------------------
# PART 2: Define functions for range generation and simulations
# ------------------------------------------------------------------------------

def generate_std_dev_analysis_routine_parameters(SIZES, MAX_SIMULATION_DEVIATION, THEORETICAL_CRITICAL_PARAMETERS):	
	
	# Initialize empty list
	STD_DEV_ANALYSIS_ROUTINE_PARAMETERS = []

	# Define left and right
	beta_c = THEORETICAL_CRITICAL_PARAMETERS["beta_c"]
	left = beta_c - MAX_SIMULATION_DEVIATION
	right = beta_c + MAX_SIMULATION_DEVIATION

	for L in SIZES:
		STD_DEV_ANALYSIS_ROUTINE_PARAMETERS.append((L, left, right))
		
	return STD_DEV_ANALYSIS_ROUTINE_PARAMETERS

def run_std_dev_analysis_simulations(TOPOLOGY, SAMPLING_PARAMETERS, STD_DEV_ANALYSIS_SAMPLING_PARAMETERS, STD_DEV_ANALYSIS_ROUTINE_PARAMETERS):
	
	Path(repo.working_tree_dir + f"/simulations/data/{TOPOLOGY}/std-dev-analysis/").mkdir(exist_ok=True) # exist_ok prevents errors if the folder
	
	# Import sampling parameters
	therm_N = SAMPLING_PARAMETERS["thermalization_samples"]
	N = SAMPLING_PARAMETERS["montecarlo_samples"]
	number_of_beta = STD_DEV_ANALYSIS_SAMPLING_PARAMETERS["number_of_plots"]
	
	# Run Markov-Chain Monte Carlo simulations
	for L, beta_min, beta_max in STD_DEV_ANALYSIS_ROUTINE_PARAMETERS:
		Path(repo.working_tree_dir + f"/simulations/data/{TOPOLOGY}/std-dev-analysis/L={L}").mkdir(exist_ok=True) # exist_ok prevents errors if the folder already exists

		i = 0
		for beta in np.linspace(beta_min, beta_max, number_of_beta):
            
            i = i+1
            print(f"\nSimulation ({i}/{number_of_beta}) for this lattice size")
		
			julia_script_filepath = repo.working_tree_dir + "/simulations/src/ising2D_cluster.jl"
			output_data_filepath = repo.working_tree_dir + f"/simulations/data/{TOPOLOGY}/std-dev-analysis/L={L}/beta={beta}.txt"
			with open(output_data_filepath, "w") as output_data_file:
				# Write the header line to the file with the current date and time
				output_data_file.write(f"# e, |m|, m2, m4 [calculated {datetime.now()} on topology: {TOPOLOGY}]\n")
			
			shell_command = "julia " + julia_script_filepath + f" {TOPOLOGY} {beta} {L} {therm_N} {N} " + output_data_filepath	
			os.system(f"{shell_command}")
			
	print("Done! Check the files at " + repo.working_tree_dir + f"/simulations/data/" + TOPOLOGY + "/std-dev-analysis/")
	return None

# ------------------------------------------------------------------------------
# PART 3: Run simulations
# ------------------------------------------------------------------------------
		
if __name__ == "__main__":

	std_dev_analysis_routine_parameters_filepath = repo.working_tree_dir + "/setup/std-dev-analysis/std_dev_analysis_routine_parameters_" + TOPOLOGY + ".txt"
	Path(repo.working_tree_dir + f"/setup/std-dev-analysis/").mkdir(exist_ok=True) # exist_ok prevents errors if the folder
	
	# Print on file
	with open(std_dev_analysis_routine_parameters_filepath,"w") as std_dev_analysis_routine_parameters_file:
		std_dev_analysis_routine_parameters_file.write(f"# L, beta_min, beta_max [calculated {datetime.now()} on topology: {TOPOLOGY}]]\n")

	STD_DEV_ANALYSIS_ROUTINE_PARAMETERS = generate_std_dev_analysis_routine_parameters(SIZES, MAX_SIMULATION_DEVIATION, THEORETICAL_CRITICAL_PARAMETERS)

	print(f"Std-dev analysis routine parameters on topology: {TOPOLOGY}")
	for i in range(len(SIZES)):
		
		# Print on file
		with open(std_dev_analysis_routine_parameters_filepath,"a") as std_dev_analysis_routine_parameters_file:
			std_dev_analysis_routine_parameters_file.write(f"{STD_DEV_ANALYSIS_ROUTINE_PARAMETERS[i][0]}, {STD_DEV_ANALYSIS_ROUTINE_PARAMETERS[i][1]}, {STD_DEV_ANALYSIS_ROUTINE_PARAMETERS[i][2]}\n")

		# Print on terminal
		interval_width = round(STD_DEV_ANALYSIS_ROUTINE_PARAMETERS[i][2] - STD_DEV_ANALYSIS_ROUTINE_PARAMETERS[i][1],3)
		print(f"{STD_DEV_ANALYSIS_ROUTINE_PARAMETERS[i]} [Interval width: {interval_width}]\n")
		
	# Comment here for debugging
	print("Running std-dev simulations...\n")
	run_std_dev_analysis_simulations(TOPOLOGY, SAMPLING_PARAMETERS, STD_DEV_ANALYSIS_SAMPLING_PARAMETERS, STD_DEV_ANALYSIS_ROUTINE_PARAMETERS)

	# Uncomment here for debugging # TODO Check, works
	# print(STD_DEV_ANALYSIS_ROUTINE_PARAMETERS)
