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
from fit_magnetic_susceptibility import quadratic_fit_function, quadratic_fit_chi, power_law_fit_function, pseudocritical_beta_fit, chi_max_fit
from plot_magnetic_susceptibility import plot_chi, plot_fss_chi

# ------------------------------------------------------------------------------
# PART 2: Define functions for critical parameters extraction
# ------------------------------------------------------------------------------

def extract_critical_parameters(TOPOLOGY, SIZES, MAX_FIT_DEVIATION, THEORETICAL_CRITICAL_PARAMETERS):

	print(f"\nStarting extraction of critical parameters on {TOPOLOGY} lattice\n")
	
	quadratic_fit_results_filepath = repo.working_tree_dir + f"/analysis/magnetic-susceptibility/results/{TOPOLOGY}/quadratic_fit_results.txt"
	with open(quadratic_fit_results_filepath,"w") as quadratic_fit_results_file:
		quadratic_fit_results_file.write(f"# L, beta_pc, e_beta_pc, chi_max, e_chi_max, chi2/ndof [calculated {datetime.now()} on topology: {TOPOLOGY}]\n")
		
	fitted_critical_parameters_filepath = repo.working_tree_dir + f"/analysis/magnetic-susceptibility/results/{TOPOLOGY}/fitted_critical_parameters.txt"
	with open(fitted_critical_parameters_filepath,"w") as fitted_critical_parameters_file:
		fitted_critical_parameters_file.write(f"# beta_c, e_beta_c, nu, e_nu, gamma/nu, e_(gamma/nu) [calculated {datetime.now()} on topology: {TOPOLOGY}]\n")
	
	quadratic_fit_results_matrix = []
	
	for L in SIZES:

		# Comment the next block if you want to manually import fit_range 
		# using customized ranges rather than the generated ones
		FIT_RANGES = generate_fit_ranges(TOPOLOGY, SIZES, MAX_FIT_DEVIATION)
		fit_range = FIT_RANGES[f"{L}"]
		
		# Uncomment the next block if you want to manually import fit_range
		# using customized ranges rather than the generated ones
		'''
		fit_ranges_filepath = repo.working_tree_dir + "/analysis/magnetic-susceptibility/fit_ranges_" + TOPOLOGY + ".txt"
		sizes, left, right = np.loadtxt(fit_ranges_filepath, delimiter=',', unpack=True)
		i = np.where(sizes==L)
		fit_range = (float(left[i]), float(right[i]))
		'''
		
		beta_pc, e_beta_pc, chi_max, e_chi_max, chi2ndof = quadratic_fit_chi(TOPOLOGY, L, fit_range)
		print(f"Quadratic fit near maximum completed. Moving on t power-law fits...\n")
		
		with open(quadratic_fit_results_filepath,"a") as quadratic_fit_results_file:
			quadratic_fit_results_file.write(f"{L}, {beta_pc}, {e_beta_pc}, {chi_max}, {e_chi_max}, {chi2ndof}\n")
			
	beta_c, e_beta_c, nu, e_nu, chi2ndof = pseudocritical_beta_fit(TOPOLOGY, SIZES, THEORETICAL_CRITICAL_PARAMETERS)
	print(f"Pseudocritical temperature power-law fit completed.")
	
	exponent, e_exponent, chi2ndof = chi_max_fit(TOPOLOGY, SIZES, THEORETICAL_CRITICAL_PARAMETERS)
	print(f"Magnetization susceptibility maxima fit completed.")
	
	gamma = exponent * nu
	e_gamma = np.sqrt((exponent * e_nu)**2 + (nu * e_exponent)**2) # TODO Not independent errors
	
	with open(fitted_critical_parameters_filepath,"a") as fitted_critical_parameters_file:
		fitted_critical_parameters_file.write(f"{beta_c}, {e_beta_c}, {nu}, {e_nu}, {exponent}, {e_exponent} \n")
		
	return None
		
# ------------------------------------------------------------------------------
# PART 3: Extraction
# ------------------------------------------------------------------------------

if __name__ == "__main__":

	extract_critical_parameters(TOPOLOGY, SIZES, MAX_FIT_DEVIATION, THEORETICAL_CRITICAL_PARAMETERS)
