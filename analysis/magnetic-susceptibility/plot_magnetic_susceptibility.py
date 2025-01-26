#!usr/bin/python3

# ------------------------------------------------------------------------------
# PART 1: Import modules and variables
# ------------------------------------------------------------------------------

import os
from pathlib import Path 		# Manage file paths
from datetime import datetime

# Get repo root
import git
repo = git.Repo('.', search_parent_directories=True)

# Import parameters
import sys
sys.path.append(repo.working_tree_dir + "/setup/")

import numpy as np
from matplotlib import pyplot as plt
plt.style.use("science")

from setup import TOPOLOGY, SIZES, THEORETICAL_CRITICAL_PARAMETERS

# ------------------------------------------------------------------------------
# PART 2: Define functions for plotting
# ------------------------------------------------------------------------------

def plot_chi(TOPOLOGY, SIZES, ax):
	
	'''
	Plot the magnetic susceptibility as a function of beta.
	Data are retrieved from the /analysis/data/ folder.
	Input:
	SIZES: (imported from setup)
	ax: figure-axis object to plot over, passed here in order to use this function elsewhere
	Output:
	ax, modified
	'''
	
	ax.set_prop_cycle('color', plt.cm.tab10.colors)

	for L in SIZES:
	
		# Load data
		input_data_filepath = repo.working_tree_dir + f"/analysis/data/{TOPOLOGY}/L={L}.txt"
		beta, _, _, chi, e_chi = np.loadtxt(input_data_filepath, delimiter=",", unpack=True)
		
		# Plot chi vs beta
		plotted_point = ax.errorbar(beta, chi, yerr=e_chi, fmt=".", label=fr"$L$ = {L}")
		ax.plot(beta, chi, "-", alpha=0.5, color=plotted_point[0].get_color()) # lines connecting the points
		
	ax.set_title(fr"$\chi'$ - {TOPOLOGY} lattice")
	ax.legend(loc="upper left")
	ax.set_xlabel(r"$\beta$")
	ax.set_ylabel(r"Magnetic susceptibility $\chi'$")
    
    # No return to act over ax

def plot_fss_chi(TOPOLOGY, SIZES, ax, beta_c, nu, exponent, e_exponent):

	'''
	Plot the finite size scaling of the magnetic susceptibility. (i.e. collapse plot).
	The data are retrieved from the /analysis/data folder.
	Input:
		SIZES: (imported from setup)
		ax: figure-axis object to plot over, passed here in order to use this function elsewhere
		beta_c: float (from theory or pseudocritical_fit)
		nu: float (from theory or pseudocritical_fit)
		gamma: float (from theory or chimax_fit)
	Output:
		ax, modified
	'''

	ax.set_prop_cycle('color', plt.cm.tab10.colors)
	    
	for L in SIZES:
    
		# Load data
		input_data_filepath = repo.working_tree_dir + f"/analysis/data/{TOPOLOGY}/L={L}.txt"
		beta, _, _, chi, e_chi = np.loadtxt(input_data_filepath, delimiter=",", unpack=True)

		# Define scaling variable and function to be plotted
		x = (beta - beta_c) * L**(1/nu)
		y = chi / L**(exponent)
		e_y = np.sqrt( ( L**(-exponent) * e_chi )**2 + ((exponent * L**(-exponent-1)) * chi * e_exponent)**2 )

		# Plot FSS function
		plotted_point = ax.errorbar(x, y, yerr=e_y, fmt=".", label=fr"$L$ = {L}")
		ax.plot(x, y, "-", alpha=0.4, color=plotted_point[0].get_color())

	ax.set_title(fr"FSS $\chi'$ - {TOPOLOGY} lattice")
	ax.legend(loc="upper right")
	ax.set_xlabel(r"$(\beta - \beta_c) L^{1/\nu}$")
	ax.set_ylabel(r"$\chi' / L^{\gamma/\nu}$")
    
    # No return to act over ax

# ------------------------------------------------------------------------------
# PART 3: Perform plots
# ------------------------------------------------------------------------------

if __name__ == "__main__":
    # Read user mode
	error_message = "No option specified \n\
Use --plot as a call option to plot data of magnetic suscpetibility obtained from processing \n\
Use --plot-fss-th as a call option to plot finite-size scaled data of magnetic susceptibility using theoretical critial parameters\n\
Use --plot-fss-fit as a call option to plot finite-size scaled data of magnetic susceptibility using critial parameters obtained from fit analysis."

	if len(sys.argv) != 2:
		raise ValueError(error_message)
	else:
		Path(repo.working_tree_dir + f"/analysis/magnetic-susceptibility/results/{TOPOLOGY}").mkdir(exist_ok=True)
		user_mode = sys.argv[1]
		if user_mode == "--plot":	
			
			# Initialize plot
			print("\nPlotting magnetic susceptibility...\n")
			
			fig, ax = plt.subplots(1, 1, figsize=(6,4.5))
			plot_chi(TOPOLOGY, SIZES, ax)
			
			plt.savefig(repo.working_tree_dir + f"/analysis/magnetic-susceptibility/results/{TOPOLOGY}/magnetic_susceptibility_plot.pdf")
			
			plt.show() # TODO Comment
			plt.close()
			
			print("Done! Check the pdf file at " + repo.working_tree_dir + f"/analysis/magnetic-susceptibility/results/{TOPOLOGY}/\n")
			
		elif user_mode == "--plot-fss-th":
			
			beta_c = THEORETICAL_CRITICAL_PARAMETERS["beta_c"]
			nu = THEORETICAL_CRITICAL_PARAMETERS["nu"]
			gamma = THEORETICAL_CRITICAL_PARAMETERS["gamma"]
			
			# Initialize plot
			print("\nPlotting FSS magnetic susceptibility using theoretical critical parameters...\n")
			
			fig, ax = plt.subplots(1, 1, figsize=(6,4.5))
			plot_fss_chi(TOPOLOGY, SIZES, ax, beta_c, nu, gamma/nu, 0.0)
			
			plt.savefig(repo.working_tree_dir + f"/analysis/magnetic-susceptibility/results/{TOPOLOGY}/magnetic_susceptibility_plot_FSS_th.pdf")
			
			plt.show()
			plt.close()
			
			print("Done! Check the pdf file at " + repo.working_tree_dir + f"/analysis/magnetic-susceptibility/results/{TOPOLOGY}/\n")
		
		elif user_mode == "--plot-fss-fit":
		
			fitted_critical_parameters_filepath = repo.working_tree_dir + f"/analysis/magnetic-susceptibility/results/{TOPOLOGY}/fitted_critical_parameters.txt"
			beta_c, _, nu, _, exponent, e_exponent = np.loadtxt(fitted_critical_parameters_filepath, delimiter=',', unpack=True)
			
			# Initialize plot
			print("\nPlotting FSS magnetic susceptibility using fitted critical parameters...\n")
			
			fig, ax = plt.subplots(1, 1, figsize=(6,4.5))
			ax.set_prop_cycle('color', plt.cm.tab10.colors)
			plot_fss_chi(TOPOLOGY, SIZES, ax, beta_c, nu, exponent, e_exponent)
			
			plt.savefig(repo.working_tree_dir + f"/analysis/magnetic-susceptibility/results/{TOPOLOGY}/magnetic_susceptibility_plot_FSS_fit.pdf")
			
			plt.show()
			plt.close()
			
			print("Done! Check the pdf file at " + repo.working_tree_dir + "/analysis/magnetic-susceptibility/results/{TOPOLOGY}/\n")
			
		else:
			raise ValueError(error_message)
