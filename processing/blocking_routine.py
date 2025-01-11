#!/usr/bin/python3

import numpy as np
import os

# Get repo root
import git
cwd = os.getcwd()
repo = git.Repo('.', search_parent_directories=True)

# Import parameters
import sys
sys.path.append(repo.working_tree_dir + "/simulations/")
from simulations_routine import parameters, number_of_beta
# from pathlib import Path # to manage file paths

# Python modules for analysis
import matplotlib.pyplot as plt
from datetime import datetime

# Define the block lengths (CHANGE!)
block_lengths = "\"1 2 3 4\""
optimal_block_length = 100

file_path_out = repo.working_tree_dir + "/processing/data/blocking_std_dev.txt" # file where data computed by blocking.jl will be stored
with open(file_path_out, "w") as file:
    # Write the header line to the file with the current date and time
    file.write(f"# L, beta, k, sigma_e, sigma_m, sigma_m2, sigma_m4 [calculated {datetime.now()}]\n")

# Run the blocking algorithm in Julia for each set of data and block length.
# The data are stored in the file /processing/data/blocking_std_dev.txt
for L, beta_min, beta_max in parameters:
    Path(f"./data/L={L}").mkdir(exist_ok=True) # exist_ok prevents errors if the folder already exists
    for beta in np.linspace(beta_min, beta_max, number_of_beta):
        shell_command = "julia " + repo.working_tree_dir + fr"/processing/src/blocking.jl 0 {L} {beta} {block_lengths}"	
        os.system(f"{shell_command}")
        
# Get back
os.chdir(cwd)

--------------------------------------------------------------------------------

# Read and plot the data to determine optimal block length
# Identical to file_path_out !!
std_dev_file = repo.working_tree_dir + "/processing/data/blocking_std_dev.txt"

L_data, beta_data, k_data, sigma_e_data, sigma_m_data, sigma_m2_data, sigma_m4_data = np.loadtxt(std_dev_file, delimiter=",", unpack=True)

def plot_std(observable):
    # Find the optimal block length. The function makes a plot with one subplot for each value
    # of L. In the subplot, the standard deviation is plotted as a function of block length, 
    # for each value of beta.

    # Choose observable
    if observable == "magnetization":
        observable = sigma_m_data
    elif observable == "magnetization2":
        observable = sigma_m2_data
    elif observable == "magnetization4":
        observable = sigma_m4_data
    elif observable == "energy":
        observable = sigma_e_data
    else:
        raise ValueError("Invalid observable. Choose 'magnetization', 'magnetization2', 'magnetization4', 'energy'.")

    plt.style.use("science")
    number_of_subplots = np.unique(L_data).size

    fig, ax = plt.subplots(number_of_subplots, 1, figsize=(8, 4*number_of_subplots), sharex=True)

    for i, L in enumerate(np.unique(L_data)): # i runs from 0 to number_of_subplots-1.
                                              # L runs over the unique values of L
        ax[i].text(0.5, 0.9, f"L = {int(L)}", transform=ax[i].transAxes, 
                   fontsize=12, horizontalalignment='center')
        ax[i].set_ylabel(r"Standard deviation $\sigma_m$")

        for beta in np.unique(beta_data):
            indices_with_L_and_beta = (L_data == L) & (beta_data == beta)
            plotted_point = ax[i].plot(k_data[indices_with_L_and_beta], observable[indices_with_L_and_beta],  
                       ".", label=fr"$\beta$ = {beta:.3f}")
            ax[i].plot(k_data[indices_with_L_and_beta], observable[indices_with_L_and_beta],  
                       "-",alpha=0.5, color=plotted_point[0].get_color()) # lines connecting the points
            ax[i].legend(loc="upper right")

    ax[-1].set_xlabel(r"Block length $k$")
    plt.savefig(repo.working_tree_dir + "/processing/data/blocking_std_dev_plot.pdf")

plot_std("magnetization2")
