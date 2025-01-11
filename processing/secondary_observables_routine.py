#!/usr/bin/python3

import numpy as np
import os
# from pathlib import Path # to manage file paths

# Get repo root
import git
cwd = os.getcwd()
repo = git.Repo('.', search_parent_directories=True)

# Import parameters
import sys
sys.path.append(repo.working_tree_dir + "/simulations/")
from simulations_routine import parameters, number_of_beta

# Run Julia resampling routine
fake_samples = 100 # Change

for L, beta_min, beta_max in parameters:
    for beta in np.linspace(beta_min, beta_max, number_of_beta):
        shell_command = "julia " + repo.working_tree_dir + fr"/processing/src/resampling.jl {L} {beta} {fake_samples}"	
        os.system(f"{shell_command}")

# Get back
os.chdir(cwd)
