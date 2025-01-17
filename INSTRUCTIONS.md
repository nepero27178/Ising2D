# A little BlaBla

The classical Ising model on 2D lattices is a good toy model to study a lot of 
statistical physics and thermodynamics. Our main goal here is to study the phase
transition occurring around the critical temperature $\beta_c$. Here all 
thermodynamic quantities undergo a second order phase transition.

Importantly, the precise value of the critical temperature depends on lattice
topology - which means that said value is different for the regular square,
triangular and hexagonal lattice. In this project we simulate all of these
options. Other quantities we try to estimate are critical exponents $\nu$ and
$\gamma$. These (and four more we are not extracting) only depend on the
broken physical symmetry when the system undergoes the phase transition. Thus
even varying the lattice topology the final results must be compatible. You will
find all details in our final relation, `Ising2D-Relation.pdf`, also available
in this repo, and in the linked bibliographyc sources.

# How do I use this repo?

This repo is organized in two layers:
1. The "bottom layer" is the one written in Julia. Scripts are contained in the
`/src/` floders, and are used to numerically simulate the lattice itself. These
scripts deal with all the numerics of the Markov-Chain Monte Carlo lattice
simulation, creating the matrix storing the lattice, updating it through a
"clustering" strategy and finally sampling it. Moreover, all the heavy-load
processing of blocking and resampling is also handled by Julia, due to the
rapidly increasing size of data. In a realistic case-of-use you will not need 
ever to run a Julia script from command line.
2. The "top layer" is the one written in Python. Apart from the `/analysis/`s and
`/setup/` folder, where everything is handled by Python rather than Julia due to
the easiness of calculations involved, in the other steps of the process Python
acts as a sort of manager: it handles all the data flow, runs the Julia scripts
in order according to the user instructions. In a realistic use-case, you will
need to call the Python scripts in a specific order we will explain in a minute,
but only setting up all details about the whole analysis (sizes to explore,
blocking and resampling parameters, and so on) at the beginning, and then
letting the entire machinery work until final results are provided.


Now, a little schematics on how all of this works.

## 1. Setup

1. The first file to worry about is `Ising2D/setup/setup.py`: by editing it
(PLEASE DO NOT CHANGE THE VARIABLES NAMES) directly from text editor, you can
set up all parameters of the simulation (like: lattice sizes to explore, lattice
topology, and so on).
2. You may have noticed some other Python scripts in the same folders: these you
can use to run a little sanity check on the parameters you have set. In
particular, `Ising2D/setup/theoretical_estimations.py` performs a very rough calculations of
where the maxima of magnetic susceptibility are expected to be for each lattice
size, with the possibility of plotting the results to let you have an idea of
how the simulations ranges should be.
3. The third file is `Ising2D/setup/routine_parameters.py` and its `.txt` analog.
By running this file you will compute for each lattice size a simulation range
centered on the expected theoretical $\beta_{pc}(L)$ (pseudocritical beta, where
the susceptibility maximum _is_) and extended left and right for an amount
decreasing as $L^{-\nu}$ (as predicted from theory). Once you have ran this
script, these parameters are all stored in `Ising2D/setup/routine_parameters.txt`
and called pretty much everywhere. *If you desire, you can change by hand the
simulated ranges* by overwriting the `.txt` file. If you do so, be careful and
run the simulations in the `--read-ranges` option (see below).

## 2. Simulations

[GO ON...]
