#!/usr/bin/python3

import numpy as np
import os
from pathlib import Path # to manage file paths

# Get repo root
import git
cwd = os.getcwd()
repo = git.Repo('.', search_parent_directories=True)

therm_N = 10000
N = 50000
number_of_beta = 10

parameters = [
#  (L, beta_min, beta_max)
    (10, 0.25, 0.29),  
    (20, 0.25, 0.29),
#    (30, 0.25, 0.29),
#	 (50, 0.25, 0.29)
]


# parameters = [(30, 0.25, 0.29)]
 
def main():  
	for L, beta_min, beta_max in parameters:
		Path(repo.working_tree_dir + f"/simulations/data/L={L}").mkdir(exist_ok=True) # exist_ok prevents errors if the folder already exists

		for beta in np.linspace(beta_min, beta_max, number_of_beta):
			shell_command = "julia " + repo.working_tree_dir + f"/simulations/src/ising2D_cluster.jl {beta} {L} {therm_N} {N} ./data/L={L}/beta={beta}.txt"	
			os.system(f"{shell_command}")
			
	os.chdir(cwd)

if __name__ == "__main__":
    main()
