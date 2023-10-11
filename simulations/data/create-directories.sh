#!/usr/bin/bash
# This code creates directories for data with different number of spins. 

for i in {2..10}; do
    N=$((2**$i))
    mkdir $N
done
