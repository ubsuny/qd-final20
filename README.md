# Eigenstates of a Quantum Dot

### Bisection method (SchrodingerSolve_user.py)

This code contains a class for solving the Schrodinger equation with the finite square well potential. The inputs for this class are:
- x: x-axis values
- b: A point outside the square well (b>>L)
- L: Size of the finite square well (-L,L)
- Vo: Potential outside the well
- en: Vector of energies used to look for eigenstates
- psi: Wavefunction values and its derivative (psi and psi')
- psi0: Initial conditions of the wavefunction

The solver class uses a bisection method for root finding and *scipy.integrate.odeint()* for solving the Schrodinger equation.

### Scipy Functions (QD_Eigenstates.ipynb)

This code contains a similar class to the previous for solving the Schrodinger equation with the finite square well potential.  The inputs for this class are the same as the previos method.  The difference in this class is the use of *scipy.optimize.brentq()* for root finding.
