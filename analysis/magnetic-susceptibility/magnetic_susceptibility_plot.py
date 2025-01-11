#!usr/bin/python3

import os
import numpy as np
from matplotlib import pyplot as plt
plt.style.use("science")

# Get repo root
import git
cwd = os.getcwd()
repo = git.Repo('.', search_parent_directories=True)

# Import parameters
import sys
sys.path.append(repo.working_tree_dir + "/simulations/")
from simulations_routine import parameters

# Set quantities estimated by the fits (exponent = gamma/nu)
EXPONENT, NU, BETA_C = 1.75, 1.0, 0.27 # CHANGE!

def plot_chi():
    """
    Plot the magnetic susceptibility as a function of beta.
    The data are retrieved from the /analysis/data folder.
    Input:
        None
    Output:
        plot saved in /analysis/magnetic-susceptibility/susceptibility.pdf
    """

    # Initialize plot
    fig, ax = plt.subplots(1, 1, figsize=(4,3) )

    for set in parameters: # i.e. for each value of L
        # Load data
        L = set[0]
        data_filepath = repo.working_tree_dir + f"/analysis/data/L={L}.txt"
        beta, _, _, chi, e_chi = np.loadtxt(data_filepath, delimiter=",", unpack=True)

        # Plot chi vs beta
        plotted_point = ax.errorbar(beta, chi, yerr=e_chi, fmt=".",
                             label=fr"$L$ = {L}")
        ax.plot(beta, chi, "-", alpha=0.5, 
             color=plotted_point[0].get_color()) # lines connecting the points

    ax.legend(loc="upper right")
    ax.set_xlabel(r"$\beta$")
    ax.set_ylabel(r"Magnetic susceptibility $\chi'$")

    plt.savefig(repo.working_tree_dir + "/analysis/magnetic-susceptibility/susceptibility.pdf")

def plot_fss_chi(exponent, nu, beta_c):
    """"
    Plot the finite size scaling of the magnetic susceptibility. (i.e. collapse plot).
    The data are retrieved from the /analysis/data folder.
    Input:
        exponent: float (gamma/nu, from chimax_fit)
        nu: float (from pseudocritical_fit)
        beta_c: float (from pseudocritical_fit)
    Output:
        plot saved in /analysis/magnetic-susceptibility/susceptibility_FSS.pdf
    """

    # Initialize plot
    plt.style.use("science")
    fig, ax = plt.subplots(1, 1, figsize=(4,3) )

    for set in parameters: # i.e. for each value of L
        # Load data
        L = set[0]
        data_filepath = repo.working_tree_dir + f"/analysis/data/L={L}.txt"
        beta, _, _, chi, e_chi = np.loadtxt(data_filepath, delimiter=",", unpack=True)

        # Define scaling variable and function to be plotted
        x_variable = (beta - beta_c) * L**(1/nu)
        y_variable = chi / L**exponent
        # TODO aggiungi errori

        # Plot FSS function
        plotted_point = ax.errorbar(x_variable, y_variable, yerr=e_chi, fmt=".", label=fr"$L$ = {L}")
        ax.plot(x_variable, y_variable, "-", alpha=0.5, color=plotted_point[0].get_color())

    ax.legend(loc="upper right")
    ax.set_xlabel(r"$(\beta - \beta_c) L^{1/\nu}$")
    ax.set_ylabel(r"$\chi' / L^{\gamma/\nu}$")

    plt.savefig(repo.working_tree_dir + "/analysis/magnetic-susceptibility/susceptibility_FSS.pdf")

def main():
    plot_chi()
    plot_fss_chi(EXPONENT, NU, BETA_C)

if __name__ == "__main__":
    main()