#!/usr/bin/bash
# This code creates directories for data with different number of spins. 

ThermSteps=1000
MonteCarloSteps=1000

for i in {2}; do
    N=$((2**$i))
    mkdir "N=$N"
    for j in {1..10}; do
    	beta=$(($j/10))
    	echo $beta
    	echo $(julia ../code/ising2D_metro.jl $beta $N $ThermSteps\
    	       $MonteCarloSteps "./N=$N/beta=$beta.txt")
    done
done
