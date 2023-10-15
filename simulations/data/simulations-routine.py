#!/usr/bin/python3 

import numpy as np
from matplotlib import pyplot as plt
import os

ThermSteps=1000
MonteCarloSteps=1000

for i in [2,3]:
    N=2**i
    os.mkdir(f"N={N}")
    for j in np.linspace(1,10,11):
    	Beta=j/10
    	Command = f"julia ../code/ising2D_metro.jl {Beta} {N} {ThermSteps} {MonteCarloSteps} ./N={N}/beta={Beta}.txt"	
    	os.system(f"{Command}")
