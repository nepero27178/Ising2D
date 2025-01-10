#!/usr/bin/python3

import numpy as np
import os
import sys
sys.path.append('../simulations/')
from simulations_routine import parameters, number_of_beta
from pathlib import Path # to manage file paths

block_lengths = "\"60 70 80 90 100 110 120\""

for L, beta_min, beta_max in parameters:
    Path(f"./data/L={L}").mkdir(exist_ok=True) # exist_ok prevents errors if the folder already exists
    for beta in np.linspace(beta_min, beta_max, number_of_beta):
        shell_command = fr"julia ./src/blocking.jl 0 {L} {beta} {block_lengths}"	
        os.system(f"{shell_command}")
