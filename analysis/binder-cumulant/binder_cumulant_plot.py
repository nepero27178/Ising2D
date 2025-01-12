#!usr/bin/python3

import os
import numpy as np
from matplotlib import pyplot as plt
# import scienceplots
# plt.style.use("science")

# Get repo root
import git
cwd = os.getcwd()
repo = git.Repo('.', search_parent_directories=True)

# Import parameters
import sys
sys.path.append(repo.working_tree_dir + "/simulations/")
from simulations_routine import parameters

def plot_U():
    """
    Plot the Binder's cumulant as a function of beta.
    The data are retrieved from the /analysis/data folder.
    Input:
        None
    Output:
        plot saved in /analysis/binder-cumulant/binder_cumulant.pdf
    """

    # Initialize plot
    fig, ax = plt.subplots(1, 1,figsize=(4,3))

    for set in parameters: # i.e. for each value of L
        # Load data
        L = set[0]
        data_filepath = repo.working_tree_dir + f"/analysis/data/L={L}.txt"
        beta, U, e_U, _, _ = np.loadtxt(data_filepath, delimiter=",", unpack=True)

        # Plot chi vs beta
        plotted_point = ax.errorbar(beta, U, yerr=e_U, fmt=".",
                             label=fr"$L$ = {L}")
        ax.plot(beta, U, "-", alpha=0.5, 
             color=plotted_point[0].get_color()) # lines connecting the points

    ax.legend(loc="upper right")
    ax.set_xlabel(r"$\beta$")
    ax.set_ylabel(r"Binder's cumulant $U$")

    plt.savefig(repo.working_tree_dir + "/analysis/binder-cumulant/binder_cumulant.pdf")

def plot_fss_U(nu, beta_c):
    """"
    Plot the finite size scaling of the Binder's cumulant (i.e. collapse plot).
    The data are retrieved from the /analysis/data folder.
    Input:
        nu: float (from pseudocritical_fit)
        beta_c: float (from pseudocritical_fit)
    Output:
        plot saved in /analysis/magnetic-susceptibility/susceptibility_FSS.pdf
    """

    # Initialize plot
    fig, ax = plt.subplots(1, 1, figsize=(4,3) )

    for set in parameters: # i.e. for each value of L
        # Load data
        L = set[0]
        data_filepath = repo.working_tree_dir + f"/analysis/data/L={L}.txt"
        beta, U, e_U, _, _ = np.loadtxt(data_filepath, delimiter=",", unpack=True)

        # Define scaling variable and function to be plotted
        x_variable = (beta - beta_c) * L**(1/nu)
        y_variable = U
        # TODO aggiungi errori

        # Plot FSS function
        plotted_point = ax.errorbar(x_variable, y_variable, yerr=e_U, fmt=".", label=fr"$L$ = {L}")
        ax.plot(x_variable, y_variable, "-", alpha=0.5, color=plotted_point[0].get_color())

    ax.legend(loc="upper right")
    ax.set_xlabel(r"$(\beta - \beta_c) L^{1/\nu}$")
    ax.set_ylabel(r"$U$")

    plt.savefig(repo.working_tree_dir + "/analysis/binder-cumulant/binder_cumulant_FSS.pdf")

def main():
	PARAMETERS = {
		# Theoretical values. TODO Update with values from pseudocritical_fit
		"nu" : 1,
		"beta_c" : np.log(3)/4
	}

	plot_U()
	plot_fss_U(PARAMETERS['nu'], PARAMETERS['beta_c'])

	return None

if __name__ == "__main__":
    main()
