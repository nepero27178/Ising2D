#!usr/bin/python3

import os
import numpy as np
from matplotlib import pyplot as plt

# Get repo root
import git
cwd = os.getcwd()
repo = git.Repo('.', search_parent_directories=True)

# Import parameters
import sys
sys.path.append(repo.working_tree_dir + "/simulations/")
from simulations_routine import parameters

# Set L

for set in parameters:
	L = set[0]

data_filepath = repo.working_tree_dir + f"/analysis/data/L={L}.txt"
beta, _, _, chi, e_chi = np.loadtxt(data_filepath)
