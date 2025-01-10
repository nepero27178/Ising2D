In this repo you will find all the files needed to simulate
the 2D Ising model on a triangular lattice. The work was related to the 
university course "Numerical Methods for Physics", held in 2023 by Professor
Claudio Bonati at the University of Pisa.

This work was done by me and a friends of mine, Marco Pompili ([@mrc-pop](https://github.com/mrc-pop)).

Three languages are used in this repo: Julia (used for simulations and low-level
numerical processing), Python (for scientific analysis of the simulated data)
and LaTeX (for the relation we wrote).

# Outline of the project
The system being simulated is a triangular lattice of $N^2$ classical spins. 
Each spin $S_{x,y}$ can assume the values $\pm 1$. The lattice shape is
rhomboidal with periodic boundary conditions, thus the lattice is encoded into a
$L \times L$ square matrix with enhanced connectivity. All details about the
model are explained inside our relation, `Ising2D-Triangular.pdf`.

This project consists of three parts in sequence: **simulation**, **processing**
and **analysis** - each of which is assigned a folder, inside of which you may
find some Python routine files used to control all the numerical machinery. In 
each you can also find a `/src/` folder and a `/data/` folder, respectively for
Julia scripts and data storage.

1. In the **simulation** part we've simulated sampling the system using a 
Markov Chain-Monte Carlo approach through the Metropolis-Hastings algorithm, in
its "clustered" form. This way the system updates its status starting from an
ordered situation (perfectly ferromagnetic), for which all the spins are 
pointing "up" (meaning, their value is $+1$); the simulation is repeated for
various choiches of $N$, spanning the temperature $\beta$. The output of each
simulation is a sampling of energy and magnetization every $N^2$ steps of the
Metropolis algorithm.

2. In the **processing** part all the datas collected are used to extract the
magnetization and therefore the magnetic susceptibility. To reduce correlation 
between data points, a "binning/blocking" scheme is applied. Lastly, some
secondary observables are useful in the last part: here we've applied the
"Boostrap" technique to get estimate statistical errors.

3. The **analysis** part is oriented to the numerical estimation of the critical
parameters in the second-order phase transition the Ising model undergoes at low
temperature. We're interested in four parameters regulating how some physical
quantities diverge at the critical temperature.

To be continued...

## References

To be continued...

# How to use this repo

To be continued...
