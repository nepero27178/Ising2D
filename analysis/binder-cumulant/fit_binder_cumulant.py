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
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
plt.style.use("science")

from setup import TOPOLOGY, SIZES, MAX_SIMULATION_DEVIATION, THEORETICAL_CRITICAL_PARAMETERS

# ------------------------------------------------------------------------------
# PART 2: Generate fit ranges
# ------------------------------------------------------------------------------

def generate_fit_range(L, NUMLEFT, NUMRIGHT):
	""" 
	Generate fit ranges for the linear fit of the binder cumulant.
	We choose NUMLEFT data points left of beta_c and NUMRIGHT right of beta_c.
	"""

	beta_c_th = THEORETICAL_CRITICAL_PARAMETERS["beta_c"]

	# Load data
	input_data_filepath = repo.working_tree_dir + f"/analysis/data/{TOPOLOGY}/L={L}.txt"
	beta, U, e_U, _, _ = np.loadtxt(input_data_filepath, delimiter=",", unpack=True)
	delta_beta = beta[1] - beta[0]

	# Find the beta_min corresponding to the NUMLEFT-th data point left of beta_c
	# and the beta_max corresponding to the NUMRIGHT-th data point right of beta_c
	beta_min = beta[beta >= beta_c_th - NUMLEFT * delta_beta][0]
	beta_max = beta[beta <= beta_c_th + NUMRIGHT * delta_beta][-1]

	return (beta_min, beta_max)

# ------------------------------------------------------------------------------
# PART 3: Fit the data 
# ------------------------------------------------------------------------------

# Step 1: linear fit

def linear_fit_function(x, a, b):
	return a * x + b

def linear_fit_U(TOPOLOGY, L, beta_min, beta_max):
	"""
	Linear fit of the Binder's cumulant U as a function of beta, 
	in the interval [beta_min, beta_max].
	"""

	# Load data
	input_data_filepath = repo.working_tree_dir + f"/analysis/data/{TOPOLOGY}/L={L}.txt"
	beta, U, e_U, _, _ = np.loadtxt(input_data_filepath, delimiter=",", unpack=True)

	# Restrict data 
	beta_fit = beta[(beta >= beta_min) & (beta <= beta_max)]
	U_fit = U[(beta >= beta_min) & (beta <= beta_max)]
	e_U_fit = e_U[(beta >= beta_min) & (beta <= beta_max)]

	# Fit
	popt, pcov = curve_fit(linear_fit_function, beta_fit, U_fit, 
						p0 = (3,-100) , sigma=e_U_fit, absolute_sigma=False)

	A, offset = popt
	e_A, e_offset = np.sqrt(np.diag(pcov))

	# Calculate chi2
	chi2 = np.sum( ((U_fit - linear_fit_function(beta_fit, *popt)) / e_U_fit)**2 )
	ndof = len(U_fit) - len(popt)
	chi2ndof = chi2 / ndof

	# PLOT
	plotted_point = ax.errorbar(beta_fit, U_fit, yerr=e_U_fit, fmt=".", label=fr"$L$ = {L}") # ($\chi^2/$ndof$ ={chi2ndof:.1f}$)")
	ax.plot(beta_fit, linear_fit_function(beta_fit, *popt), "-", alpha=0.5, color=plotted_point[0].get_color()) # best-fit lines
	ins_ax.plot(beta_fit, linear_fit_function(beta_fit, *popt), "-", alpha=0.5, color=plotted_point[0].get_color()) # best-fit lines

	return A, e_A, offset, e_offset, chi2ndof

# Step 2: power-law fit

def power_law_fit_function(x, a, b):
	return a * np.power(x, b)

def power_law_fit_U(TOPOLOGY):
	"""
	Power-law fit of the Binder's cumulant's slope as a function of L.
	"""

	# Read results for A(L) from linear fit, for all sizes
	input_data_filepath = repo.working_tree_dir + f"/analysis/binder-cumulant/results/{TOPOLOGY}/linear_fit_results.txt"
	L, A, e_A, _, _, _  = np.loadtxt(input_data_filepath, delimiter=",", unpack=True)

	A = np.abs(A) # Making A positive just for convenience

	# Fit
	popt, pcov = curve_fit(power_law_fit_function, L, A, 
						p0=(-1.0,1.0), sigma=e_A, absolute_sigma=True)
	
	abs_Uprimestar, invnu = popt
	e_abs_Uprimestar, e_invnu = np.sqrt(np.diag(pcov))

	nu = 1 / invnu
	e_nu = e_invnu / invnu**2 

	# Calculate chi2
	chi2 = np.sum( ((A - power_law_fit_function(L, *popt)) / e_A)**2 )
	ndof = len(L) - len(popt)
	print(f"chi2 = {chi2}")
	chi2ndof = chi2 / ndof

	# PLOT
	plotted_point = ax.errorbar(L, A, yerr=e_A, fmt=".", color="tab:blue", label="Estimated slopes at intersection")
	ax.plot(L, power_law_fit_function(L, *popt), "-", color="tab:cyan", alpha=0.5, label=rf"Best-fit function ($\chi^2$/ndof $={chi2ndof:.2f}$)") # best-fit lines

	return abs_Uprimestar, e_abs_Uprimestar, nu, e_nu, chi2ndof


def estimate_intersection(TOPOLOGY):
	
	# Read results for A(L) from linear fit, for all sizes
	input_data_filepath = repo.working_tree_dir + f"/analysis/binder-cumulant/results/{TOPOLOGY}/linear_fit_results.txt"
	L, A, e_A, offset, e_offset, _  = np.loadtxt(input_data_filepath, delimiter=",", unpack=True)

	# Define matrix and vector for LLS
	M = np.ones([len(L),2], float)
	M[:,1] = - A
	b = offset

	# the intersection point minimizes ||Mx - b||
	intersection, e_intersection, _, _ = np.linalg.lstsq(M, b, rcond=None)

	return intersection, e_intersection

# ------------------------------------------------------------------------------
# PART 4: Main run
# ------------------------------------------------------------------------------

if __name__ == "__main__":

	error_message = "No option specified \n\
Use --fit-linear as a call option for the linear fit near the intersection (Step 1)\n\
Use --fit-power-law as a call option to fit A(L) vs L (Step 2) \n\
Use --intersection as a call option to estimate intesection (Step 3) \n\
"
	if len(sys.argv) != 2:
		raise ValueError(error_message)
	else:
		user_mode = sys.argv[1]
			
		if user_mode == "--fit-linear":	

			user_input = input("Insert NUMLEFT, NUMRIGHT: ")

			NUMLEFT, NUMRIGHT = tuple(int(x) for x in user_input.split(","))

			# Initialize file with results
			linear_fit_results_filepath = repo.working_tree_dir + f"/analysis/binder-cumulant/results/{TOPOLOGY}/linear_fit_results.txt"
			with open(linear_fit_results_filepath,"w") as linear_fit_results_file:
				linear_fit_results_file.write(f"# L, A, e_A, offset, e_offset, chi2/ndof [calculated {datetime.now()} on topology: {TOPOLOGY}, NUMLEFT={NUMLEFT}, NUMRIGHT={NUMRIGHT}] \n")
		
			fig, ax = plt.subplots(figsize=(6, 4.5))

			ax.set_prop_cycle('color', plt.cm.tab10.colors)
			ax.set_title(r"Linear fits on $U$ near the intersection")
			ax.set_xlabel(r"$\beta$")
			ax.set_ylabel(r"$U$")
			
			ins_ax = ax.inset_axes([.65, .55, .3, .4])  # [x, y, width, height] w.r.t. ax

			# Fit and write on file
			for L in SIZES:
				beta_min, beta_max = generate_fit_range(L, NUMLEFT, NUMRIGHT)
				A, e_A, offset, e_offset, chi2ndof = linear_fit_U(TOPOLOGY, L, beta_min, beta_max)
				with open(linear_fit_results_filepath,"a") as linear_fit_results_file:
					linear_fit_results_file.write(f"{L}, {A}, {e_A}, {offset}, {e_offset}, {chi2ndof}\n")
				print(f"Linear fit for L={L}, topology={TOPOLOGY} completed.\n")
		
			intersection, e_intersection = estimate_intersection(TOPOLOGY)
			#ax.errorbar(intersection[1], intersection[0], yerr=e_intersection, xerr=e_intersection, color="black", label="Intersection")
			ax.plot(intersection[1], intersection[0], "*", markersize=5, color="black", label="Intersection")
			ins_ax.plot(intersection[1], intersection[0], "*",  markersize=5, color="black", label="Intersection")


			#axins = ax.inset_axes([0.2752, 1.175, 7e-4, 0.022], xlim=(0.2743, 0.275), ylim=(1.155, 1.175), transform=ax.transData)
			if TOPOLOGY == "triangular":
				ins_ax.set_xlim([0.27455, 0.27475])
				ins_ax.set_ylim([1.163, 1.172])
			elif TOPOLOGY == "hexagonal":
				ins_ax.set_xlim([0.6582, 0.6590])
				ins_ax.set_ylim([1.199, 1.213])
			elif TOPOLOGY == "square":
				ins_ax.set_xlim([0.44, 0.45])
				ins_ax.set_ylim([1.17, 1.18])
			
			ins_ax.set_xticklabels([])
			ins_ax.set_yticklabels([])
			ax.indicate_inset_zoom(ins_ax, ls="-", edgecolor="gray", alpha=0.5)

			ax.legend(loc="lower left")

			plt.savefig(repo.working_tree_dir + f"/analysis/binder-cumulant/results/{TOPOLOGY}/binder_linear_fit_intersection.pdf")
			plt.show()


			print(f"Fits completed! The results have been written on file.")


		elif user_mode == "--fit-power-law":
			print(f"Fitting power law...")

			fig, ax = plt.subplots(figsize=(6,4.5))

			ax.set_prop_cycle('color', plt.cm.tab10.colors)
			ax.set_xlabel(r"$L$")
			ax.set_ylabel(r"$|A(L)|$")
			ax.set_title(r"Power-law fit of Binder cumulant's slope")
			abs_Uprimestar, e_abs_Uprimestar, nu, e_nu, chi2ndof = power_law_fit_U(TOPOLOGY)
			print(f"Power-law fit completed.\n\
Fitted parameters:\n\
|U_2'(0)| = {abs_Uprimestar} +/- {e_abs_Uprimestar}\n\
nu = {nu} +/- {e_nu}\n\
chi^2/ndof = {chi2ndof}")
			
			ax.legend()
			plt.savefig(repo.working_tree_dir + f"/analysis/binder-cumulant/results/{TOPOLOGY}/binder_power_law_fit.pdf")
			plt.show()
		
		elif user_mode == "--intersection":
			print(f"Reading data from file, and estimating intersection...")
			intersection, e_intersection = estimate_intersection(TOPOLOGY)
			print(f"Done! Estimated intersection: (x,y)=({intersection[1]},{intersection[0]}), error={e_intersection[0]}")

		else:
			raise ValueError(error_message)