#!/usr/bin/python3

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

# ------------------------------------------------------------------------------
# PART 2: Check duality relation
# ------------------------------------------------------------------------------

fitted_critical_parameters_triangular_filepath = repo.working_tree_dir + "/analysis/magnetic-susceptibility/results/triangular/fitted_critical_parameters.txt"
bt, e_bt, _, _, _, _ = np.loadtxt(fitted_critical_parameters_triangular_filepath, delimiter=',', unpack=True)
# bt, e_bt = (0.27465, 5e-5) # From Binder

fitted_critical_parameters_hexagonal_filepath = repo.working_tree_dir + "/analysis/magnetic-susceptibility/results/hexagonal/fitted_critical_parameters.txt"
bh, e_bh, _, _, _, _ = np.loadtxt(fitted_critical_parameters_hexagonal_filepath, delimiter=',', unpack=True)
# bh, e_bh = (0.6585, 1e-4) # From Binder

constraint = np.sinh(2*bt) * np.sinh(2*bh)
e_constraint = 2 * np.sqrt( (np.cosh(2*bt) * np.sinh(2*bh) * e_bt)**2 + (np.sinh(2*bt) * np.cosh(2*bh) * e_bh)**2 )

print(f"Duality relation constraint (product of hypersines): {constraint} +/- {e_constraint}")
