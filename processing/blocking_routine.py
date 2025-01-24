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

# Import parameters
import sys
sys.path.append(repo.working_tree_dir + "/setup/")
from setup import TOPOLOGY, SIZES, SAMPLING_PARAMETERS, STD_DEV_ANALYSIS_SAMPLING_PARAMETERS

# Read which simulations have been performed
routine_parameters_filepath = repo.working_tree_dir + "/setup/routine_parameters_" + TOPOLOGY + ".txt"
_, left, right = np.loadtxt(routine_parameters_filepath, delimiter=',', unpack=True)
ROUTINE_PARAMETERS = []

for i in range(len(SIZES)):
	ROUTINE_PARAMETERS.append((SIZES[i], left[i], right[i]))
	
del left, right

# Read which simulations for std-dev analysis have been performed
std_dev_analysis_routine_parameters_filepath = repo.working_tree_dir + "/setup/std-dev-analysis/std_dev_analysis_routine_parameters_" + TOPOLOGY + ".txt"
_, left, right = np.loadtxt(std_dev_analysis_routine_parameters_filepath, delimiter=',', unpack=True)
STD_DEV_ANALYSIS_ROUTINE_PARAMETERS = []

for i in range(len(SIZES)):
	STD_DEV_ANALYSIS_ROUTINE_PARAMETERS.append((SIZES[i], left[i], right[i]))
	
del left, right
		
# Python modules for analysis
import matplotlib.pyplot as plt
plt.style.use("science")
from datetime import datetime
		
# ------------------------------------------------------------------------------
# PART 2: Define functions for blocking
# ------------------------------------------------------------------------------

def try_blocking_lengths(TOPOLOGY, STD_DEV_ANALYSIS_ROUTINE_PARAMETERS, STD_DEV_ANALYSIS_SAMPLING_PARAMETERS):

	'''
	try_blocking_lengths runs the julia blocking.jl script over each blocking
	length specified in \"Ising2D/setup/setup.pt\". These lenghts are imported,
	converted to string, and passed to the julia script as input arguments. The
	result of the double for cycle is a single file, at output_data_filepath,
	where for each L, beta and block length k standard deviations of different
	observables are saved.
	'''
	
	print(f"\nTrying different blocking schemes on the {TOPOLOGY} lattice...\n")

	# Import samping parameters
	number_of_beta = STD_DEV_ANALYSIS_SAMPLING_PARAMETERS["number_of_plots"]
	block_trial_lengths = STD_DEV_ANALYSIS_SAMPLING_PARAMETERS["block_trial_lengths"]
	string = ""
	
	for length in block_trial_lengths:
		string = string + f"{length} "
		
	block_trial_lengths = "\"" + string[:-1] + "\""
	
	# TODO Find a clever way to indicate filepaths
	
	output_data_filepath = repo.working_tree_dir + f"/processing/std-dev-analysis/{TOPOLOGY}/blocking_std_dev.txt" # file where data computed by blocking.jl will be stored
	Path(repo.working_tree_dir + "/processing/std-dev-analysis/").mkdir(exist_ok=True)
	Path(repo.working_tree_dir + f"/processing/std-dev-analysis/{TOPOLOGY}").mkdir(exist_ok=True)
	
	with open(output_data_filepath, "w") as output_data_file:
		# Write the header line to the file with the current date and time
		output_data_file.write(f"# L, beta, k, sigma_e, sigma_|m|, sigma_m2, sigma_m4 [calculated {datetime.now()} on topology: {TOPOLOGY}]\n")

	# Run the blocking algorithm in Julia for each set of data and block length.
	# The data are stored in the file /processing/data/trial_blocking_std_dev.txt
	for L, beta_min, beta_max in STD_DEV_ANALYSIS_ROUTINE_PARAMETERS:
	
		print(f"Currently working on block size {L}")		
		julia_script_filepath = repo.working_tree_dir + "/processing/src/blocking.jl"				# julia script filepath
		
		for beta in np.linspace(beta_min, beta_max, number_of_beta):
		
			input_data_filepath = repo.working_tree_dir + f"/simulations/data/{TOPOLOGY}/std-dev-analysis/L={L}/beta={beta}.txt" # Use raw data
			shell_command = "julia " + julia_script_filepath + f" {L} {beta} {block_trial_lengths} " + input_data_filepath + " " + output_data_filepath + " --try"
			os.system(f"{shell_command}")
			
	print("\nDone! Check the file at " + output_data_filepath + "\n")
	return None
	
def plot_std_dev(TOPOLOGY, STD_DEV_ANALYSIS_ROUTINE_PARAMETERS, STD_DEV_ANALYSIS_SAMPLING_PARAMETERS):

	'''
	plot_std_dev plots the data calculated by try_blocking_lengths by directly
	reading the saved file.
	'''
	
	# Import samping parameters
	number_of_beta = STD_DEV_ANALYSIS_SAMPLING_PARAMETERS["number_of_plots"]
	beta_min = STD_DEV_ANALYSIS_ROUTINE_PARAMETERS[0][1]	# All equal
	beta_max = STD_DEV_ANALYSIS_ROUTINE_PARAMETERS[0][2]	# All equal
	
	std_dev_filepath = repo.working_tree_dir + f"/processing/std-dev-analysis/{TOPOLOGY}/blocking_std_dev.txt" # file where data computed by blocking.jl will be stored
	sizes, betas, lengths, sigma_e, sigma_m, sigma_m2, sigma_m4 = np.loadtxt(std_dev_filepath, delimiter=',', unpack=True)
	sigma = [sigma_e, sigma_m, sigma_m2, sigma_m4]
	observables_names = [r"$\sigma_e$", r"$\sigma_{|m|}$", r"$\sigma_{m^2}$", r"$\sigma_{m^4}$"]
	
	for beta in np.linspace(beta_min, beta_max, number_of_beta):
		
		fig, ax = plt.subplots(len(sigma), 1, figsize=(8, 4*len(sigma)), sharex=True)
		
		for L in SIZES:
		
			indices = (sizes == L) & (betas == beta)
		
			for i in range(len(sigma)):

				std_dev = sigma[i]
				ax[i].set_ylabel(r"Standard deviation " + observables_names[i])

				plotted_point = ax[i].plot(lengths[indices], std_dev[indices], ".", label=fr"$L$ = {L}")	# Points
				ax[i].plot(lengths[indices], std_dev[indices], "-",alpha=0.5, color=plotted_point[0].get_color()) 	# Lines connecting the points
				ax[i].set_title(fr"std:dev from {TOPOLOGY} lattice at $\beta$ = {beta:.2f}")
				ax[i].legend(loc="upper right")
		
		plt.savefig(repo.working_tree_dir + f"/processing/std-dev-analysis/{TOPOLOGY}/blocking_std_dev_plot_beta={beta}.pdf")
		plt.close()
	
	print("\nDone! Check the plots at " + repo.working_tree_dir + f"/processing/std-dev-analysis/{TOPOLOGY}/\n")
	return None
	
def block_data(TOPOLOGY, ROUTINE_PARAMETERS, SAMPLING_PARAMETERS):

	'''
	blocking_data runs the julia blocking.jl script just once for the optimal 
	block length in \"Ising2D/setup/setup.pt\". The result of the double for 
	cycle is a system of files specular to that of \"Ising2D/simulations/data\",
	at output_data_filepath, where for each L, and beta sampled values are
	blocked and saved for the specified optimal block length.
	'''

	# Import sampling parameters
	number_of_beta = SAMPLING_PARAMETERS["number_of_beta"]
	block_optimal_length = SAMPLING_PARAMETERS["block_optimal_length"]
	Path(repo.working_tree_dir + f"/processing/data/{TOPOLOGY}/blocklength={block_optimal_length}").mkdir(exist_ok=True)
	
	print(f"\nBlocking data from the {TOPOLOGY} lattice using optimal block length: {block_optimal_length}\n")

	# Run the blocking algorithm in Julia for each set of data and block length.
	# The data are stored in the file /processing/data/trial_blocking_std_dev.txt
	for L, beta_min, beta_max in ROUTINE_PARAMETERS:
		
		for beta in np.linspace(beta_min, beta_max, number_of_beta):
		
			input_data_filepath = repo.working_tree_dir + f"/simulations/data/{TOPOLOGY}/L={L}/beta={beta}.txt"										# Use raw data
			output_data_filepath = repo.working_tree_dir + f"/processing/data/{TOPOLOGY}/blocklength={block_optimal_length}/L={L}/beta={beta}.txt" 	# Path to file where data computed by blocking.jl will be stored
			
			Path(repo.working_tree_dir + f"/processing/data/{TOPOLOGY}/blocklength={block_optimal_length}/L={L}").mkdir(exist_ok=True)
			with open(output_data_filepath, "w") as output_data_file:
				# Write the header line to the file with the current date and time
				output_data_file.write(f"# e, |m|, m2, m4 [calculated {datetime.now()} on topology: {TOPOLOGY}]\n")
				
			julia_script_filepath = repo.working_tree_dir + "/processing/src/blocking.jl"			# Run julia script
			shell_command = "julia " + julia_script_filepath + f" {L} {beta} {block_optimal_length} " + input_data_filepath + " " + output_data_filepath + " --use-optimal"
			os.system(f"{shell_command}")
	
	print("\nDone! Check the files at " + repo.working_tree_dir + f"/processing/data/{TOPOLOGY}/blocklength={block_optimal_length}/\n")
	return None

# ------------------------------------------------------------------------------
# PART 3: Run optimal blocking OR trial blocking
# ------------------------------------------------------------------------------

if __name__ == "__main__":

	# Read user mode
	error_message = "No option specified \n\
Use --use-optimal as a call option to use the previously computed optimal block length \n\
Use --try as a call option to try different block lengths \n\
Use --plot as a call option to plot previously calculated standard deviations of observables \n\
Use --try-plot as a call option to try different block lengths and plot standard deviations of observables \n\
In order to change the used optimal block length, or to set trial lengths, modify their values in \"Ising2D/setup/setup.py\"."

	if len(sys.argv) != 2:
		raise ValueError(error_message)
	else:
		user_mode = sys.argv[1]
		if user_mode == "--try":
			
			# Comment here for debugging
			try_blocking_lengths(TOPOLOGY, STD_DEV_ANALYSIS_ROUTINE_PARAMETERS, STD_DEV_ANALYSIS_SAMPLING_PARAMETERS)
			
			# Uncomment here for debugging
			# print(ROUTINE_PARAMETERS)
			
		elif user_mode == "--plot":
			
			# Comment here for debugging
			plot_std_dev(TOPOLOGY, STD_DEV_ANALYSIS_ROUTINE_PARAMETERS, STD_DEV_ANALYSIS_SAMPLING_PARAMETERS)
			
			# Uncomment here for debugging
			# print(ROUTINE_PARAMETERS)
		
		elif user_mode == "--try-plot":
			
			# Comment here for debugging
			try_blocking_lengths(TOPOLOGY, STD_DEV_ANALYSIS_ROUTINE_PARAMETERS, STD_DEV_ANALYSIS_SAMPLING_PARAMETERS)
			plot_std_dev(TOPOLOGY, STD_DEV_ANALYSIS_ROUTINE_PARAMETERS, STD_DEV_ANALYSIS_SAMPLING_PARAMETERS)
			
			# Uncomment here for debugging
			# print(ROUTINE_PARAMETERS)
			
		elif user_mode == "--use-optimal":

			# Comment here for debugging
			Path(repo.working_tree_dir + f"/processing/data/{TOPOLOGY}/").mkdir(exist_ok=True)
			block_data(TOPOLOGY, ROUTINE_PARAMETERS, SAMPLING_PARAMETERS)
			
			# Uncomment here for debugging
			# print(ROUTINE_PARAMETERS)
			
		else:
			raise ValueError(error_message)
