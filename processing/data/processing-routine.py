#!/usr/bin/python3

import numpy as np
from matplotlib import pyplot as plt
import os

for i in [1,2,3,4]:
    L=10*i
    os.mkdir(f"L={L}")
    for j in np.linspace(4,5,11):
    	Beta=format(j/10,'.4f')
    	Command = f"julia ../code/blocking.jl {Beta} {L} ../../simulations/data/L={L}/beta={Beta}.txt ./L={L}/beta={Beta}.txt"
    	os.system(f"{Command}")