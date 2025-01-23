In this repo you will find all the files needed to simulate
the 2D Ising model on a triangular lattice. The work was related to the 
university course "Numerical Methods for Physics", held in 2023 by Professor
Claudio Bonati at the University of Pisa.

This work was done by me and a friend of mine, Marco Pompili ([@mrc-pop](https://github.com/mrc-pop)).

Three languages are used in this repo: Julia (used for simulations and low-level
numerical processing), Python (for scientific analysis of the simulated data)
and LaTeX (for the final report we wrote).

# Outline of the project
The system being simulated is a triangular or hexagonal lattice of $L^2$ classical 
spins. Each spin $s_{\mathbf{x},\mathbf{y}}$ can take the values $\pm 1$. The lattice 
shape is rhomboidal with periodic boundary conditions, thus the lattice is encoded 
into a $L \times L$ square matrix with increased or decreased connectivity (the 
coordination number, i.e. the number of nearest neighbours for each site, is equal 
to $6$ for the triangular lattice and $3$ for the hexagonal lattice). All details 
about the model are explained inside our relation, `Ising2D_hexagonal_triangular.pdf`.

This project consists of four sequential parts: **setup**, **simulations**, 
**processing** and **analysis**. Each of them is assigned a folder, inside of
which you may find some Python routine files used to control all the numerical
machinery. In each you can also find a `/src/` folder and a `/data/` folder, 
respectively for Julia scripts and data storage.

0. The **setup** step is self-explanatory: here we choose which sizes and temperatures
to use for our subsequent simulations, as well as other important choices for the 
analysis, for instance, the blocking length. Additionally, we store the theoretical
estimations for the critical exponents and other quantities related to the phase
transition of the model.

1. In the **simulations** part we've simulated sampling the system using a 
Markov Chain-Monte Carlo approach through the Metropolis-Hastings algorithm, with
its "cluster" variant (Wolff's algorithm). This way the system updates its status 
starting from an ordered situation (perfectly ferromagnetic), for which all the 
spins are pointing "up" (meaning, their value is $+1$); the simulation is repeated 
for various choiches of $L$, spanning the temperature $\beta$. The output of each
simulation is a sampling of energy, magnetization (plus its square and fourth power).

2. In the **processing** part all the collected data are used to extract the
magnetization and therefore the magnetic susceptibility. To reduce correlation 
between data points, a "binning/blocking" scheme is applied. Lastly, some
secondary observables are useful in the last part: here we've applied the
"Boostrap" technique to estimate their statistical errors appropriately.

3. The **analysis** part is oriented to the numerical estimation of the critical
parameters in the second-order phase transition the Ising model undergoes at low
temperature. We're interested in four parameters regulating how some physical
quantities diverge at the critical temperature. In particular, we study the 
(modified) magnetic susceptivity $\chi'$ and the Binder cumulant $U$, in order
to get an estimate for the critical exponents and the critical temperature.

To be continued...

## References

To be continued...

# How to use this repo

To be continued...
