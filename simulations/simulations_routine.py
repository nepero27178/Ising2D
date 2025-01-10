#!/usr/bin/python3 

import numpy as np
import os
from pathlib import Path # to manage file paths

# Go to right directory
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

therm_N = 10000
N = 500
number_of_beta = 4

parameters = [
#  (L, beta_min, beta_max)
    (10, 0.25, 0.29),  
    (15, 0.25, 0.29),
#   (30, 0.25, 0.29),
]


# parameters = [(30, 0.25, 0.29)]

 
def main():  
	for L, beta_min, beta_max in parameters:
		Path(f"./data/L={L}").mkdir(exist_ok=True) # exist_ok prevents errors if the folder already exists

		for beta in np.linspace(beta_min, beta_max, number_of_beta):
			shell_command = f"julia ./src/ising2D_cluster.jl {beta} {L} {therm_N} {N} ./data/L={L}/beta={beta}.txt"	
			os.system(f"{shell_command}")

if __name__ == "__main__":
    main()
