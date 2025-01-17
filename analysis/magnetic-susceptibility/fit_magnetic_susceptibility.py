#!usr/bin/python3

# ------------------------------------------------------------------------------
# PART 1: Import modules and variables
# ------------------------------------------------------------------------------

import os
from datetime import datetime

# Get repo root
import git
repo = git.Repo('.', search_parent_directories=True)

# Import parameters
import sys
sys.path.append(repo.working_tree_dir + "/setup/")

import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
# plt.style.use("science")

from setup import TOPOLOGY, SIZES, MAX_FIT_DEVIATION, THEORETICAL_CRITICAL_PARAMETERS
from generate_fit_ranges import generate_fit_ranges
	
# Import plot function for plottind data along with fits
# from plot_magnetic_susceptibility import plot_chi, plot_fss_chi
    
# ------------------------------------------------------------------------------
# PART 2: Define functions for fitting
# ------------------------------------------------------------------------------
    
# Quadratic fit function
def quadratic_fit_function(x, a, b, c):
	return a * (x - b)**2 + c
    
def quadratic_fit_chi(TOPOLOGY, L, fit_range):

	'''
	Quadratic fit for the magnetic susceptibility, near the maximum. 
	The point of maximum is by definition beta_pc(L), pseudocritical.
	Beta and L are considered exact parameters.
	Input:
		TOPOLOGY: square, triangular or hexagonal?
		L: lattice size
		fit_range: for given size, the selected fit range
	Output:
		best parameters, with errors, and chi^2.
	'''
    
	# Load data
	input_data_filepath = repo.working_tree_dir + f"/analysis/data/{TOPOLOGY}/L={L}.txt"
	beta, _, _, chi, e_chi = np.loadtxt(input_data_filepath, delimiter=",", unpack=True)
	
	# Restrict data
	left, right = fit_range
	
	processed_points = len(beta)
	deleted_points = 0
	deleted_indices = np.where((beta<left)|(beta>right))
	beta = np.delete(beta,deleted_indices)
	chi = np.delete(chi,deleted_indices)
	e_chi = np.delete(e_chi,deleted_indices)
	deleted_points = len(deleted_indices[0])
		
	simulation_to_fit_efficiency = (processed_points-deleted_points)/processed_points
	if simulation_to_fit_efficiency < 0.25:
		print(f"\nHey! You are discarding a lot of data here! [{TOPOLOGY} lattice, L={L}, efficiency={simulation_to_fit_efficiency:.2f}]\n")
    
	# Import theoretically estimated parameters to initialize p0
	theoretical_critical_parameters_filepath = repo.working_tree_dir + "/setup/theoretical_estimations_" + TOPOLOGY + ".txt"
	sizes, tmp_beta_pc, tmp_chi_max = np.loadtxt(theoretical_critical_parameters_filepath, delimiter=',', unpack=True)
	
	not_found = True
	for index in range(len(sizes)):
		if L==sizes[index]:
			not_found = False
			th_beta_pc = tmp_beta_pc[index]
			th_chi_max = tmp_chi_max[index]
	
	if not_found:
		AssertionError(f"Seems like the chosen size L={L} was not simulated in the first place.\n\
Please run at least one simulation and the processing procedure to perform this fit!\n")
		sys.exit()
    
	# Fit (-1 arbitrarily set due to negative concavity)
	popt, pcov = curve_fit(quadratic_fit_function, beta, chi, p0=(-1, th_beta_pc, th_chi_max), sigma=e_chi, absolute_sigma=True)

	# Extract parameters (A is discarded)
	_, beta_pc, chi_max = popt
	_, e_beta_pc, e_chi_max = np.sqrt(np.diag(pcov))
	
	# Compute chi^2
	chi2 = np.sum( ((chi - quadratic_fit_function(beta, *popt)) / e_chi)**2 )
	ndof = len(beta)-3
	chi2ndof = chi2/ndof
	
	return beta_pc, e_beta_pc, chi_max, e_chi_max, chi2ndof

# Power law fit function
def power_law_fit_function(x, a, b, c):
	
	return a + b * np.power(x, c)

def pseudocritical_beta_fit(TOPOLOGY, SIZES, THEORETICAL_CRITICAL_PARAMETERS):
    
	'''
	Fit the pseudocritical beta as a function of L.
	Input:
		TOPOLOGY: square, triangular or hexagonal?
		SIZES: to check if all available data are being fitted
		THEORETICAL_CRITICAL_PARAMETERS: to import fit initializers
	Output:
		best parameters, with errors, and chi^2.
	'''

	# Load data
	input_data_filepath = repo.working_tree_dir + f"/analysis/magnetic-susceptibility/results/{TOPOLOGY}/quadratic_fit_results.txt"
	L, beta_pc, e_beta_pc, _, _, _ = np.loadtxt(input_data_filepath, delimiter=",", unpack=True)
	
	if len(L) < len(SIZES):
	
		# TODO Improve to L != SIZES
		print(f"\nWarning: I have detected you are trying to run a fit over {len(L)+1} data points, while you have set a longer array SIZES in \"Ising2D/setup/setup.py\". You may be fitting over an undernumbered array of points. Consider running the quadratic fit over the remaining lengths in order to increment the fit quality.")
		user_proceed = input("Shall I proceed?  (y/n) ")
		
		if user_proceed == "n":
			sys.exit()
		elif user_proceed == "y":
			print("Fitting undernumbered data...\n")
		else:
			raise ValueError("Invalid input. Please enter y (yes) or n (no).")

	th_beta_c = THEORETICAL_CRITICAL_PARAMETERS["beta_c"]	# Theoretical critical temperature
	th_nu = THEORETICAL_CRITICAL_PARAMETERS["nu"]			# Theoretical value of exponent nu
	th_x0 = 1 												# x0 set arbitrarily to 1
	
	th_exponent = -1/th_nu

	# Fit
	popt, pcov = curve_fit(power_law_fit_function, L, beta_pc, p0=(th_beta_c, th_x0, th_exponent),
		        sigma=e_beta_pc, absolute_sigma=True)

	# Extract parameters (x0 is discarded)
	beta_c, _, invnu = popt
	e_beta_c, _, e_invnu = np.sqrt(np.diag(pcov))

	nu = 1 / invnu
	e_nu = e_invnu / invnu**2 # error propagation

	# Calculate chi^2
	chi2 = np.sum( ((beta_pc - power_law_fit_function(L, *popt)) / e_beta_pc)**2 )
	ndof = len(beta)-3
	chi2ndof = chi2/ndof
	
	return beta_c, e_beta_c, nu, e_nu, chi2ndof
	
def chi_max_fit(TOPOLOGY, SIZES, THEORETICAL_CRITICAL_PARAMETERS):
    
	'''
	Fit the susceptibility maxima as a function of L.
	Input:
		TOPOLOGY: square, triangular or hexagonal?
		SIZES: to check if all available data are being fitted
		THEORETICAL_CRITICAL_PARAMETERS: to import fit initializers
	Output:
		best parameters, with errors, and chi^2.
	'''

	# Load data
	input_data_filepath = repo.working_tree_dir + f"/analysis/magnetic-susceptibility/results/{TOPOLOGY}/quadratic_fit_results.txt"
	L, _, _, chi_max, e_chi_max, _ = np.loadtxt(input_data_filepath, delimiter=",", unpack=True)
	
	if len(L) < len(SIZES):
	
		# TODO Improve to L != SIZES
		print(f"\nWarning: I have detected you are trying to run a fit over {len(L)+1} data points, while you have set a longer array SIZES in \"Ising2D/setup/setup.py\". You may be fitting over an undernumbered array of points. Consider running the quadratic fit over the remaining lengths in order to increment the fit quality.")
		user_proceed = input("Shall I proceed?  (y/n) ")
		
		if user_proceed == "n":
			sys.exit()
		elif user_proceed == "y":
			print("Fitting undernumbered data...")
		else:
			raise ValueError("Invalid input. Please enter y (yes) or n (no).")

	th_nu = THEORETICAL_CRITICAL_PARAMETERS["nu"]			# Theoretical value of exponent nu
	th_gamma = THEORETICAL_CRITICAL_PARAMETERS["gamma"]		# Theoretical value of exponent gamma
	th_y0 = 1 												# y0 (prefactor) set arbitrarily to 1
	th_c0 = 0												# c0 (offset) set arbitrarily to 1
	
	th_exponent = th_gamma/th_nu

	# Fit
	popt, pcov = curve_fit(power_law_fit_function, L, beta_pc, p0=(th_c0, th_y0, th_exponent), sigma=e_chi_max, absolute_sigma=True)

	# Extract parameters (x0 is discarded)
	_, _, exponent = popt
	_, _, e_exponent = np.sqrt(np.diag(pcov))

	# Calculate chi^2
	chi2 = np.sum( ((beta_pc - power_law_fit_function(L, *popt)) / e_beta_pc)**2 )
	ndof = len(beta)-3
	chi2ndof = chi2/ndof
	
	return exponent, e_exponent, chi2ndof

# ------------------------------------------------------------------------------
# PART 3: Run calculations
# ------------------------------------------------------------------------------

if __name__ == "__main__":
	# Read user mode
	error_message = "No option specified \n\
Use --fit-chi as a call option to fit data on magnetic susceptibility (chi) via a quadratic function near the maximum \n\
Use --fit-beta_pc as a call option to fit data on magnetic susceptibility maximum position beta_pc via a power-law function \n\
Use --fit-chi_max as a call option to fit data on magnetic susceptibility maximum chi_max via a power-law function \n\
For each option except --fit-all you will be asked for a specific size L; through these options it is possible to plot the resulting fit.\n\
To change the used fit ranges modify [...] END"

	if len(sys.argv) != 2:
		raise ValueError(error_message)
	else:
		user_mode = sys.argv[1]
			
		if user_mode == "--fit-chi":
		
			try:
				L = int(input("Specify desired length for data to be fitted: "))
				print(f"Near-maximum quadratic fitting for L={L}...")
			except ValueError:
				print("Invalid input! Please, enter an integer input.")
				sys.exit()
			
			# Comment the next block if you want to manually import fit_range 
			# using customized ranges rather than the generated ones
			FIT_RANGES = generate_fit_ranges(TOPOLOGY, SIZES, MAX_FIT_DEVIATION)
			fit_range = FIT_RANGES[f"{L}"]
			
			# Uncomment the next block if you want to manually import fit_range
			# using customized ranges rather than the generated ones
			'''
			fit_ranges_filepath = repo.working_tree_dir + "/analysis/magnetic-susceptibility/fit_ranges.txt"
			sizes, left, right = np.loadtxt(fit_ranges_filepath, delimiter=',', unpack=True)
			i = np.where(sizes==L)
			fit_range = (float(left[i]), float(right[i]))
			'''
			
			beta_pc, e_beta_pc, chi_max, e_chi_max, chi2ndof = quadratic_fit_chi(TOPOLOGY, L, fit_range)
			print(f"Quadratic fit near maximum completed.\n\
Fitted parameters [L={L}, fit range={fit_range}]:\n\
beta_pc = {beta_pc} +/- {e_beta_pc}\n\
chi_max = {chi_max} +/- {e_chi_max}\n\
chi^2/ndof = {chi2ndof}")

		elif user_mode == "--fit-beta_pc":
		
			print(f"Pseudocritical temperature power-law fitting...")
			beta_c, e_beta_c, nu, e_nu, chi2ndof = pseudocritical_beta_fit(TOPOLOGY, SIZES, THEORETICAL_CRITICAL_PARAMETERS)
			print(f"Pseudocritical temperature power-law fit completed.\n\
Fitted parameters:\n\
beta_c = {beta_c} +/- {e_beta_c}\n\
nu = {nu} +/- {e_nu}\n\
chi^2/ndof = {chi2ndof}")

		elif user_mode == "--fit-chi_max":
		
			print(f"Magnetization susceptibility maxima power-law fitting...")
			exponent, e_exponent, chi2ndof = chi_max_fit(TOPOLOGY, SIZES, THEORETICAL_CRITICAL_PARAMETERS)
			print(f"Magnetization susceptibility maxima fit completed.\n\
Fitted parameter:\n\
gamma/nu = {exponent} +/- {e_exponent}\n\
chi^2/ndof = {chi2ndof}")
			
		else:
			raise ValueError(error_message)
