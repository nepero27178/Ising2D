In this repo you will find all the files needed to simulate
the 2D Ising model on a square lattice. The work was related to the university
course "Numerical Methods for Physics", held in 2023 by Professor Claudio
Bonati in the University of Pisa.

This work was done by me and some friends of mine:
- Valerio Caporioni
- Alessandro Gori (me)
- Marco Pompili
- (Jacopo Resasco)

You will probably find similar works on their repos.

The vast majority of the code is written in C and Julia, but we've also used
Python and Bash to analyze the data. You will also find a LaTeX document I
wrote as a relation on my numerical experiment.

# Outline of the project
This project consists of three main parts: **simulation**, **analysis** and
**critical phenomena**. The system being simulated is a square lattice of $N^2$ 
spins $1/2$. Every spin $S_{x,y}$ can assume the values $\pm 1$.

1. In the **simulation** part we've simulated sampling the system using a 
Markov chain approach through the Metropolis algorithm. This way the system
updates its status starting from an ordered situation, where all the spins are
pointing "up" (meaning, their value is $+1$); the simulation is repeated for
various choiches of $N$, spanning the temperature $\beta$. The output of each
simulation is a sampling of energy and magnetization every $N^2$ steps of the
Metropolis algorithm.

2. In the **analysis** part all the datas collected are used to extract the
magnetization and therefore the magnetical susceptibility. To reduce correlation 
between data points, a "binning/blocking" scheme is applied. Lastly, some
secondary observables are useful in the last part: here we've applied both the
"boostrap" and the "jack-knife" techniques.

3. The **critical phenomena** part is oriented to the numerical approximation
of the critical parameters in the second-order phase transition the Ising
model undergoes at low temperature. We're interested in four parameters regulating
how some physical quantities diverge at the critical temperature. We've followed
two different paths: the **Binder's cumulant** can provide a first approximation
for $\beta_c$, the critical temperature, and for the critical exponent $\nu$.
These are used as guesses to fit the magnetic susceptibility in a **finite**
**size scaling** framework.

## References

To be continued...

# How to use this repo

To be continued...
