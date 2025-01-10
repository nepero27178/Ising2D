#!/usr/bin/python3

import numpy as np
import os
import sys
sys.path.append('../simulations/')
from simulations_routine import parameters, number_of_beta
from pathlib import Path # to manage file paths

optimal_block_length = 100 # change

for L, beta_min, beta_max in parameters:
    Path(f"./data/L={L}").mkdir(exist_ok=True) # exist_ok prevents errors if the folder already exists
    for beta in np.linspace(beta_min, beta_max, number_of_beta):
        shell_command = f"julia ./src/blocking.jl 1 {L} {beta} {optimal_block_length}"	
        os.system(f"{shell_command}")
