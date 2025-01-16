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
# plt.style.use("science")

from setup import SIZES, MAX_FIT_DEVIATION

# ------------------------------------------------------------------------------
# PART 2: Generate fit ranges
# ------------------------------------------------------------------------------

def generate_fit_ranges(SIZES, MAX_FIT_DEVIATION):

	beta_pc_filepath = repo.working_tree_dir + "/setup/theoretical_estimations.txt"
	_, beta_pc, _ = np.loadtxt(beta_pc_filepath, delimiter=',', unpack=True)

	FIT_RANGES = {}

	for i in range(len(SIZES)):

		beta_min = beta_pc[i] - MAX_FIT_DEVIATION * np.min(SIZES)/SIZES[i]
		beta_max = beta_pc[i] + MAX_FIT_DEVIATION * np.min(SIZES)/SIZES[i]
		FIT_RANGES[f"{SIZES[i]}"] = (beta_min, beta_max)
		
	return FIT_RANGES
	
# ------------------------------------------------------------------------------
# PART 3: Main run, write on file
# ------------------------------------------------------------------------------

if __name__ == "__main__":
	
	fit_ranges_filepath = repo.working_tree_dir + "/analysis/magnetic-susceptibility/fit_ranges.txt"
	with open(fit_ranges_filepath, "w") as fit_ranges_file:
		fit_ranges_file.write(f"# L, beta_min, beta_max [calculated {datetime.now()}]\n")
		
	FIT_RANGES = generate_fit_ranges(SIZES, MAX_FIT_DEVIATION)
	for L in SIZES:
		beta_min, beta_max = FIT_RANGES[f"{L}"]
		with open(fit_ranges_filepath, "a") as fit_ranges_file:
			fit_ranges_file.write(f"{L}, {beta_min}, {beta_max}\n")
