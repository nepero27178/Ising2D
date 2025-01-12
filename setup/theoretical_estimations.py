#!/usr/bin/python3
'''
This script computes the expected theoretical values for the magnetization
maxima in the 2D triangular Ising model. The dictionary PARAMETERS stores the
first-order values of the model. From the module routine_parameters.py, SIZES
are imported: there are the lattice sizes the user simulates.

Run this script with options --hide-ranges not to plot the beta ranges of the
simulations. Run this script with option --plot-ranges to plot the beta ranges
of the simulations.
'''

import os
import sys

import numpy as np
from matplotlib import pyplot as plt
# import scienceplots
# plt.style.use("science") # TODO style

import random

# Get repo root
import git
cwd = os.getcwd()
repo = git.Repo('.', search_parent_directories=True)

# Import parameters
from setup import SIZES

# Expected magnetic susceptibility parameters
PARAMETERS = {
	# Physical parameters
	"beta_c" : np.log(3)/4,
	"nu" : 1,
	"gamma" : 7/4,
	# Numerical parameters
	"x_0" : -0.1,
	"y_0" : 0.1,
	"c_0" : 0
}

# Functions

def get_theoretical_values(PARAMETERS,SIZES):
	'''
	Expected maxima of the magnetic susceptibility are stored inside the file
	theoretical_estimations.txt, and calculated via the finite size scaling
	relations of tht model.
	'''

	filepath = repo.working_tree_dir + "/setup/theoretical_estimations.txt"

	# Re-initialize file
	with open(filepath,"w") as file:
	    file.write("# L, beta_pc, chi_max\n")
	
	# Physical parameters
	beta_c = PARAMETERS['beta_c']
	nu = PARAMETERS['nu']
	gamma = PARAMETERS['gamma']

	# Numerical parameters
	x_0 = PARAMETERS['x_0'] 
	y_0 = PARAMETERS['y_0']
	c_0 = PARAMETERS['c_0']

	# Initialize theoretical values
	theoretical_values = np.zeros([len(SIZES),3])

	for i in range(len(SIZES)):
		L = SIZES[i]
		theoretical_values[i,0] = L
		theoretical_values[i,1] = beta_c + x_0 * L**(-1/nu)
		theoretical_values[i,2] = c_0 + y_0 * L**(gamma/nu)

		# Write on file
		with open(filepath,"a") as file:
		    file.write(f"{L}, {theoretical_values[i,1]}, {theoretical_values[i,2]}\n")

		'''
		# Print on terminal
		print("--------------------------------------\n")
		print(f"beta_pc[L={L}] = {round(theoretical_values[i,0],5)}\n")
		print(f"chi_max[L={L}] = {round(theoretical_values[i,1],5)}\n")
		'''
		
	return theoretical_values
	
def plot_theoretical_values(theoretical_values,ranges):
	'''
	Plot the values calculated via get_theoretical_values; together, are also
	plotted the simulations ranges stored in the imported variable RANGES
	'''
	
	beta_pc = theoretical_values[:,1]
	chi_max = theoretical_values[:,2]

	# Initialize plot
	fig, ax = plt.subplots(1, 1, figsize=(8,6))
		
	# Draw lines	
	ax.plot(beta_pc, chi_max, "-", alpha=0.5, color="gray")	
	
	for i in range(len(SIZES)):
	
		# Draw each point with the relative simulations range (define color via r)
		L = int(SIZES[i])
		r = np.linspace(0,1,len(SIZES))[i]
		
		if ranges:
			left_side = ranges[i][1]
			right_side = ranges[i][2]
			ax.plot([left_side, right_side],[chi_max[i], chi_max[i]],'-',color=(r,1-r,r))
			
		# Plot chi_max vs beta_pc
		ax.plot(beta_pc[i], chi_max[i],'.',label=fr"$L$ = {L}",color=(r,1-r,r))
	
	ax.plot(np.log(3)/4 * np.ones(len(chi_max)),
			np.linspace(0,np.max(chi_max),len(chi_max)),
			color="gray",alpha=0.5,
			linewidth=0.75)
	
	# General styling
	ax.set_xticks(ticks=[0.266, 0.270, 0.274, np.log(3)/4, 0.278],
				  labels=["0.266", "0.270", "0.274", r"$\beta_c$", "0.278"])
	ax.set_title("Theoretical values")
	ax.legend(loc="upper left")
	ax.set_xlabel(r"$\beta$")
	ax.set_ylabel(r"Magnetic susceptibility maximum $\chi_\max'$")

	# Save plot
	plt.savefig(repo.working_tree_dir + "/setup/routine_ranges_and_theoretical_maxima.pdf")

	return None

if __name__ == "__main__":
	theoretical_values = get_theoretical_values(PARAMETERS,SIZES)
	
	if len(sys.argv)!=2:
		print("Error: no option specified \n\
Use --plot-ranges as a call option to plot simulations ranges \n\
Use --hide-ranges as a call option not to plot simulations ranges \n")
		raise AssertionError()
	else:
		use_ranges = sys.argv[1]
		if use_ranges == "--plot-ranges":
			from routine_parameters import ranges
			plot_theoretical_values(theoretical_values,ranges)
		elif use_ranges == "--hide-ranges":
			ranges = False
			plot_theoretical_values(theoretical_values,ranges)
		else:
			print("Error: no option specified \n\
Use --plot-ranges as a call option to plot simulations ranges \n\
Use --hide-ranges as a call option not to plot simulations ranges \n")
			raise AssertionError()
