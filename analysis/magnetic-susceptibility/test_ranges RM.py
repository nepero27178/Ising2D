#!usr/bin/python3

import os
import numpy as np
from matplotlib import pyplot as plt
#plt.style.use("science")

# Get repo root
import git
cwd = os.getcwd()
repo = git.Repo('.', search_parent_directories=True)

# Import parameters
import sys
sys.path.append(repo.working_tree_dir + "/simulations/")
from simulations_routine import parameters

# Import parameters
sys.path.append(repo.working_tree_dir + "/setup/")
from setup import SIZES
from routine_parameters import ranges

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
	

def plot_chi_and_ranges(theoretical_values,ranges):

	beta_pc = theoretical_values[:,1]
	chi_max = theoretical_values[:,2]

	# Initialize plot
	fig, ax = plt.subplots(1, 1, figsize=(8,6))
		
	# Draw lines	
	ax.plot(beta_pc, chi_max, "-", alpha=0.5, color="gray")	
	
	for i in range(len(SIZES)):
	
		# Draw each point with the relative simulations range (define color via r)
		L = int(SIZES[i])
		
		data_filepath = repo.working_tree_dir + f"/analysis/data/L={L}.txt"
		beta, _, _, chi, e_chi = np.loadtxt(data_filepath, delimiter=",", unpack=True)
		
		# Plot chi vs beta
		plotted_point = ax.plot(beta, chi,".", label=fr"$L$ = {L}")
		ax.plot(beta, chi, "-", alpha=0.5, color=plotted_point[0].get_color()) # lines connecting the points
		
		if ranges:
			left_side = ranges[i][1]
			right_side = ranges[i][2]
			ax.plot([left_side, right_side],[chi_max[i], chi_max[i]],'-', color=plotted_point[0].get_color())
			
		# Plot chi_max vs beta_pc
		ax.plot(beta_pc[i], chi_max[i],'.',label=fr"$L$ = {L}", color=plotted_point[0].get_color())

	ax.legend(loc="upper right")
	ax.set_xlabel(r"$\beta$")
	ax.set_ylabel(r"Magnetic susceptibility $\chi'$")
	
	ax.plot(np.log(3)/4 * np.ones(len(chi_max)),
		np.linspace(0,np.max(chi_max),len(chi_max)),
		color="gray",alpha=0.5,
		linewidth=0.75)

	plt.savefig(repo.working_tree_dir + "/analysis/magnetic-susceptibility/test_ranges.pdf")

if __name__ == "__main__":
	theoretical_values = get_theoretical_values(PARAMETERS,SIZES)
	plot_chi_and_ranges(theoretical_values,ranges)
