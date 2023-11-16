#!/usr/bin/python3 

import numpy as np
from matplotlib import pyplot as plt
import os

ThermSteps=1000
MonteCarloSteps=1000

for i in [1,2,3,4]:
    N=10*i
    os.mkdir(f"N={N}")
    for j in np.linspace(4,5,11):
    	Beta=format(j/10,'.4f')
    	Command = f"julia ../code/ising2D_metro.jl {Beta} {N} {ThermSteps} {MonteCarloSteps} ./N={N}/beta={Beta}.txt"	
    	os.system(f"{Command}")
