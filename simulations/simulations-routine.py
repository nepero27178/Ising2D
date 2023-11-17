#!/usr/bin/python3 

import numpy as np
import os

ThermSteps=1000
MonteCarloSteps=1000

for i in [1,2,3,4]:
    L=10*i
    os.mkdir(f"L={L}")
    for j in np.linspace(4,5,11):
    	Beta=format(j/10,'.4f')
    	Command = f"julia ./code/ising2D_metro.jl {Beta} {L} {ThermSteps} {MonteCarloSteps} ./L={L}/beta={Beta}.txt"	
    	os.system(f"{Command}")
