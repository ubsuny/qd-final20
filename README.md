# Eigenstates of a Quantum Dot

### Finite Difference method (SchrodingerSolve_user.py)

This code contains a class for solving the Schrodinger equation with the finite square well potential. The inputs for this class are:
- x: x-axis values
- b: A point outside the square well (b>>L)
- L: Size of the finite square well (-L,L)
- Vo: Potential outside the well
- en: Vector of energies used to look for eigenstates
- m: Mass of the particle in the well.  For quantum dots, the mass of an electron is used.
- hbar: The reduced planck constant (h/2*pi)

The solver class uses a finite difference method for numerically solving the Schrodinger equation.  An analytical solution, using a root finding method, is also included to check numerical solutions.


### Scipy Functions (SchrodingerSolver_scipy.py)

This code contains a similar class to the previous for solving the Schrodinger equation with the finite square well potential.  The function *scipy.integrate.ode()* is used to numerically solve the Schrodinger equation.  Two additional inputs are required to run this class:
- psi: Wavefunction values and its derivative (psi and psi')
- psi0: Initial conditions of the wavefunction

### Jupyter Notebook (QD_Eigenstates.ipynb)

This notebook is used to compare the numerical solutions from the two solver classes to the analyitcal solution.

### Unit Test (unitTest.py)

The unit test checks to see if the wavefunctions determined by the finite difference code are normalized, i.e.,

$$ \int_{\infty}^{\infty} |\Psi^2| = 1$$
