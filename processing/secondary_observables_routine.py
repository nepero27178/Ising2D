#!/usr/bin/python3

import numpy as np
import os
from datetime import datetime
# from pathlib import Path # to manage file paths

# Get repo root
import git
cwd = os.getcwd()
repo = git.Repo('.', search_parent_directories=True)

# Import parameters
import sys
sys.path.append(repo.working_tree_dir + "/simulations/")
from simulations_routine import parameters, number_of_beta

# Number of fake samples (R) to perform
fake_samples = 100 # Change

for L, beta_min, beta_max in parameters:

    # Write the header line to the file with the current date and time
    file_path_out = repo.working_tree_dir + f"/analysis/data/L={L}.txt" # file where data computed by resampling.jl will be stored
    with open(file_path_out, "w") as file:
        file.write(f"# beta, U, e_U, chi e_chi [calculated {datetime.now()}]\n")

    # Go to processing folder
    os.chdir(repo.working_tree_dir + "/processing")

    # Run the julia resampling script
    for beta in np.linspace(beta_min, beta_max, number_of_beta):
        shell_command = "julia " + repo.working_tree_dir + fr"/processing/src/resampling.jl {L} {beta} {fake_samples}"	
        os.system(f"{shell_command}")

# Get back
os.chdir(cwd)
