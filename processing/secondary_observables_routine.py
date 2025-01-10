#!/usr/bin/python3

import numpy as np
import os
import sys
sys.path.append('../simulations/')
from simulations_routine import parameters, number_of_beta
from pathlib import Path # to manage file paths

fake_samples = 100 # Change

for L, beta_min, beta_max in parameters:
    for beta in np.linspace(beta_min, beta_max, number_of_beta):
        shell_command = fr"julia ./src/resampling.jl {L} {beta} {fake_samples}"	
        os.system(f"{shell_command}")
