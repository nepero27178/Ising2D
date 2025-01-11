#!usr/bin/python3

import os
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit

# Get repo root
import git
cwd = os.getcwd()
repo = git.Repo('.', search_parent_directories=True)

# Import parameters
import sys
sys.path.append(repo.working_tree_dir + "/simulations/")
from simulations_routine import parameters

def quadratic_fit_chi(beta, chi, e_chi):
    """
    Quadratic fit for the magnetic susceptibility, near the maximum. 
    The point of maximum is by definition beta_pc(L), pseudocritical.
    Beta and L are considered exact parameters.
    Input:
        beta: np.array
        chi: np.array
        e_chi: np.array
    Output:
        best parameters, with errors, and chi^2.
    """
    # Fit function
    def fit_func(x, a, b, c):
        return a * (x - b)**2 + c
    
    # Fit
    popt, pcov = curve_fit(fit_func, beta, chi, p0=(-1, 0.27, np.max(chi)),
                           sigma=e_chi, absolute_sigma=True)
    
    # Extract parameters (A is discarded)
    _, beta_pc, chi_max = popt
    _, e_beta_pc, e_chi_max = np.sqrt(np.diag(pcov))

    # calculate chi^2
    chi2 = np.sum( ((chi - fit_func(beta, *popt)) / e_chi)**2 )

    return beta_pc, e_beta_pc, chi_max, e_chi_max, chi2


def pseudocritical_fit(L, beta_pc, e_beta_pc):
    """
    Fit the pseudocritical beta as a function of L.
    Input:
        L: np.array
        beta_pc: np.array (estimated by quadratic_fit_chi)
        e_beta_pc: np.array (estimated by quadratic_fit_chi)
    Output:
        best parameters, with errors, and chi^2.
    """

    # Fit function
    def fit_func(x, a, b, c):
        return a + b * np.power(x, -c)

    # Fit
    popt, pcov = curve_fit(fit_func, L, beta_pc, p0=(0.27, 1, 1), # x0 set arbitrarily to 1
                        sigma=e_beta_pc, absolute_sigma=True)

    # Extract parameters (x0 is discarded)
    beta_c, _, invnu = popt
    e_beta_c, _, e_invnu = np.sqrt(np.diag(pcov))

    nu = 1 / invnu
    e_nu = e_invnu / invnu**2 # error propagation

    # Calculate chi^2
    chi2 = np.sum( ((beta_pc - fit_func(L, *popt)) / e_beta_pc)**2 )

    return beta_c, e_beta_c, nu, e_nu, chi2

def chimax_fit(L, chi_max, e_chi_max):
    """
    Fit the maximum of the magnetic susceptibility as a function of L.
    Input:
        L: np.array
        chi_max: np.array (estimated by quadratic_fit_chi)
        e_chi_max: np.array (estimated by quadratic_fit_chi)
    Output:
        best parameters, with errors, and chi^2.
    """
    # Fit function TODO valutare se serve un terzo parametro di offset
    def fit_func(x, a, b):
        return a * x**b
    
    # Fit
    popt, pcov = curve_fit(fit_func, L, chi_max, p0=(1, 7/4), # y0 set arbitrarily to 1
                           sigma=e_chi_max, absolute_sigma=True)
    
    # Extract parameters. Exponent is our estimate of gamma/nu
    y0, exponent = popt
    e_y0, e_exponent = np.sqrt(np.diag(pcov))

    # Calculate chi^2
    chi2 = np.sum( ((chi_max - fit_func(L, *popt)) / e_chi_max)**2 )

    return y0, e_y0, exponent, e_exponent, chi2

def calculate_gamma(nu, e_nu, exponent, e_exponent):
    """
    Calculate gamma from the two fit results nu and exponent=gamma/nu.
    """
    gamma = nu * exponent
    e_gamma = 0.1  ### np.sqrt((exponent * e_nu)**2 + (nu * e_exponent)**2) # TODO AMMESSO CHE SIANO INDIPENDENTI (DUBITO!) DA SISTEMARE!!
    return gamma, e_gamma

def main():

    # Define a three-column matrix to store L, beta_pc, e_beta_pc.
    # It will be needed for the pseudocritical fit.
    first_fit_results_matrix = []

    for set in parameters: # i.e. for each value of L
        # Load data
        L = set[0]
        data_filepath = repo.working_tree_dir + f"/analysis/data/L={L}.txt"
        beta, _, _, chi, e_chi = np.loadtxt(data_filepath, delimiter=",", unpack=True)

        # Quadratic fit
        beta_pc, e_beta_pc, chi_max, e_chi_max, chi2 = quadratic_fit_chi(beta, chi, e_chi)
        print(f"Quadratic fit for L={L} completed. " +
              f"beta_pc = {beta_pc:.4f} +/- {e_beta_pc:.4f}, " + 
              f"chi_max = {chi_max:.4f} +/- {e_chi_max:.4f}, " +
              f"chi^2 = {chi2:.4f}.\n")

        # Store the results
        first_fit_results_matrix.append([L, beta_pc, e_beta_pc, chi_max, e_chi_max])

    array_L, array_beta_pc, array_e_beta_pc, array_chi_max, array_e_chi_max = np.array(first_fit_results_matrix).T

    # Pseudocritical fit
    beta_c, e_beta_c, nu, e_nu, chi2 = pseudocritical_fit(array_L, array_beta_pc, array_e_beta_pc)
    print(f"Pseudocritical fit for L={L} completed. " +
            f"beta_c = {beta_c:.4f} +/- {e_beta_c:.4f}, " + 
            f"nu = {nu:.4f} +/- {e_nu:.4f}, " +
            f"chi^2 = {chi2:.4f}.\n")
    # TODO per ora non riesce a convergere... i dati sono troppo pochi

    # Chi max fit
    y0, e_y0, exponent, e_exponent, chi2 = chimax_fit(array_L, array_chi_max, array_e_chi_max)
    print(f"Chi max fit for L={L} completed. " +
            f"y0 = {y0:.4f} +/- {e_y0:.4f}, " + 
            f"exponent = {exponent:.4f} +/- {e_exponent:.4f}, " +
            f"chi^2 = {chi2:.4f}.\n")

    # Estimate gamma from the fits
    gamma, e_gamma = calculate_gamma(nu, e_nu, exponent, e_exponent)
    print(f"Estimate of gamma for L={L} is {gamma:.4f} +/- {e_gamma:.4f}.\n")

if __name__ == "__main__":
    main()