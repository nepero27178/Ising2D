The 2D classical Ising model for 1/2 spin on a square lattice is a good toy 
model to study a lot of statistical physics and thermodynamics. Our main goal 
is to study the phase transition occurring around the critical temperature
$T_c \approx 2.27 K$. Here all thermodynamic quantities undergo a second order
phase transition.
Critical exponents, typical of power-law dependencies in physical quantities
around the critcal point, can be studied in the Ising model. Here we assumed the
hamiltonian to be the simplest possible, with a ferromagnetic exchange
interaction between nearest neighbours
$$
H = -J \sum_{\langle x,y \rangle} S_x S_y
\qquad
J > 0
$$
The mathematical details regarding all that follows can be found in our document
*[insert relation]*. Here we only expose our line of reasoning to extract all
the investigated parameters.

1. **Lattice simulations**: We use a Metropolis algorithm to simulate a square
lattice of $L \times L$ classical spins. $L$ is varied from $10$ to $100$ with
steps of $10$. For each $L$, $\beta$ is varied in $21$ points around the
critical temperature (which is about %\beta_c \approx 0.44/k_B$). [to do...]
