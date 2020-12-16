# Quantum Dot eigenstates

### 1. Introduction

#### a. Quantum Dots

Quantum dots (QDs) are nanometer-scale semiconductor materials that, due to their size, display quantum confinement, in which electrons cannot escape the "dot".  This behavior allows for the "particle in a box" approximation to be used to model the QD system.  Quantum dots have been a topic of high interest in recent years due to their numerous applications, such as LEDs, single-electron transistors and solar cells.

### b. Finite Potential Square Well

The Finite Potential Square Well is a well studied solution to the Schrodinger equation for which there exists an analytical solution.  The system is time independent and therefore the Schordinger equation simplifies to its time-independent form,

$$ \frac{-\hbar^2}{2m} \frac{d^2 \Psi \left( x \right)}{dx^2} + V \left( x \right) \Psi \left( x \right) = E \Psi \left( x \right). $$

The potential energy function for the finite potential square well is,

$$ V \left( x \right) =
\begin{cases}
	0, & for -a \leq x \leq a \\
	V_0, & for \left| x \right| > a \\
\end{cases}
$$

where $V_0$ is a positive constant.  Therefore, the solution breaks into a system of three equations defined for each region.

For the regions $x < -a$ and $x > a$, the Schrodinger equation becomes,

$$ \frac{-\hbar^2}{2m} \frac{d^2 \Psi \left( x \right)}{dx^2} + V_0 \Psi \left( x \right) = E \Psi \left( x \right). $$

This can be simplified into a general, solvable ODE,

$$ \frac{d^2 \Psi \left( x \right)}{dx^2} = \alpha^2 \Psi \left( x \right), $$

where,

$$ \alpha = \frac{\sqrt{2m \left( V_0 - E \right)}}{\hbar^2}. $$

This ODE has a well known general solution,

$$ \begin{aligned}
\Psi_1 \left( x \right) &= A e^{-\alpha x} + B e^{\alpha x} \\
\Psi_3 \left( x \right) &= F e^{-\alpha x} + G e^{\alpha x}
\end{aligned}$$

where $\Psi_1 \left(x\right)$ is the solution for the region $x<-a$ and $\Psi_3 \left(x\right)$ is the solution for the region $x>a$.  For $x<-a$, $A e^{-\alpha x}$ blows up as $x \rightarrow-\infty$, and only the second term is considered.  Likewise, for $x>a$, $G e^{\alpha x}$ blows up as $x \rightarrow \infty$, so only the first term is considered.

For the region $-a \leq x \leq a$, the Schrodinger equation becomes

$$ \frac{-\hbar^2}{2m} \frac{d^2 \Psi \left( x \right)}{dx^2} = E \Psi \left( x \right). $$

This can similarly be simplified to the solvable ODE,

$$\frac{d^2 \Psi \left( x \right)}{dx^2} = \kappa^2 \Psi \left( x \right), $$

where,

$$ \alpha = \frac{\sqrt{-2m E}}{\hbar^2}. $$

The general solution for this ODE is the same as above, however for this section, a solution with trigonometric functions is considered.  This form of the general solution is then,

$$ \Psi_2 \left( x \right) = C sin\left(\kappa x\right) + D cos\left(\kappa x\right). $$

This equation can further divide the solution into two cases: symmetric (even-parity) and antisymmertic (odd-parity).  For the symmetric case, we consider only the term $D cos\left(\kappa x\right)$ and for the antisymmetric case the term $C sin\left(\kappa x\right)$.

Applying continuity conditions $\Psi_1 \left( x = -a \right) = \Psi_2 \left( x = -a \right)$, $\Psi_3 \left( x = a \right) = \Psi_2 \left( x = a \right)$, $\frac{d\Psi_1}{dx} \left( x = -a \right) = \frac{d\Psi_2}{dx} \left( x = -a \right)$, $\frac{d\Psi_3}{dx} \left( x = a \right) = \frac{d\Psi_2}{dx} \left( x = a \right)$ leads to the transendental equations,

$$ \begin{aligned}
tan\left(z\right) &= \sqrt{\left(z_0/z\right)^2 - 1} \;\;\; \left(symmetric \, case\right) \\
-cot\left(z\right) &= \sqrt{\left(z_0/z\right)^2 - 1} \;\;\; \left(symmetric \, case\right)
\end{aligned}$$

where $z_0 = \frac{a}{\hbar^2}\sqrt{2m V_0}$ and $z = \frac{a}{\hbar}\sqrt{2m\left(V_0-E\right)}$.  Thus to find analytical solutions for the eigenstates of the finite potential square well, numerical or graphical methods are used to find the roots of the transendental equations.


### 2. Computational Methods

To solve the Schrodinger equation computationally, multiple approaches can be taken.  Possibly the most straightforward approach is to use numerical root finding functions on the transendental equations from the previous section to find the energy values corresponding to each eigenstate.  Given the eigenstates and the analytical solution, wavefunctions for each state can be found.


Another approach to finding eigenstates is to find all wavefunctions in a given energy range that converge to zero outside the potential well.  In more depth, values for $\Psi\left(x\rightarrow\infty \right)$ are found for a range of energies such that $E \leq V_0$.  The root finding algorithm is then used to find those energies for which $\Psi\left(x\rightarrow\infty \right) = 0$.  These energies represent the bounded states of this system where the wavefunction does not diverge outside the potential well. Using this method is effective rather accurate but requires two numerical algorithms, root finding and integration, which can be costly depending on the system.  This method is explored in the provided code using functions from the Python library Scipy.

The final approach considered for this problem is to use a Finite Difference discretization of the Schrodinger equation to find the eigenvalues and thus the eigenstates of the system.  This method is more versatile because the system of discretized equations can be derived directly from the general Schrodinger equation and it only requires an algorithm to find eigenvalues.

#### a. Root finding

The Python library Scipy contains the function *scipy.optimize.brentq()* which is based on Brent's method for root finding.  Brent's method is a hybrid method that combines the bisection  and secant method.  Unlike the Secant method, which uses a linear interpolation curve, Brent's method uses quadratic interpolations to find the roots of the given function.

#### b. Integration

##### Backwards Differentiaion/Adams Method

The Python libray Scipy also contains the function *scipy.integrate.odeint()*, which uses the FORTRAN library odepack to solve initial value problems.  Specifically, for stiff systems, the Adams method is used and for non-stiff systems, the Backwards Differentiation method (BDM) is used.  As these methods solve initial value problems, the Schrodinger equation must be formulated into a state-space representation with given initial values for $\Psi|_{\left(x=0\right)}$ and $\frac{d\Psi}{dx}|_{\left(x=0\right)}$.  The state-space representation of the Schrodinger equation used in this code is given by,

$$\begin{aligned}
\dot \psi &= \frac{d\Psi\left(x\right)}{dx} \\
\ddot \psi &= \frac{2*m}{\hbar^2}(V(x) - E)*\Psi\left(x\right).
\end{aligned}$$

##### Finite Difference Method

The technique primarily used in this code is a finite difference method which utilizes approximations of derivatives to obtain solutions.  It is worthwhile to simplify the Schrodinger equation here for use in this method.  The wavefuntion, $\Psi \left( x \right)$, is factored out of the left hand side of the equation giving,

$$ \left[\frac{-\hbar^2}{2m} \frac{d^2}{dx^2} + V \left( x \right) \right] \Psi \left( x \right) = E \Psi \left( x \right). $$

The expression $H = \frac{-\hbar^2}{2m} \frac{d^2}{dx^2} + V \left( x \right)$ is called the *Hamiltonian operator* and allows the equation to be simplified to the form,

$$ H \Psi \left( x \right) = E \Psi \left( x \right). $$

The Hamiltonian operator is commonly expressed as the sum of the kinetic and potential energy, $H = T + V$.  Thus, $\frac{-\hbar}{2m} \frac{d^2 \Psi \left( x \right)}{dx^2}$ is equivalent to the kinetic energy.

The first derivative of the wavefunction is given in two forms, forward and backwards differences,

$$ \begin{aligned}
\Psi_{forward}'\left(x\right) &= \frac{\Psi\left(x+h\right) - \Psi\left(x\right)}{h}, \\
\Psi_{backward}'\left(x\right) &= \frac{\Psi\left(x\right) - \Psi\left(x-h\right)}{h},
\end{aligned} $$

where h is the step size.  From these expressions, the second derivative is approximated using a central difference method,

$$ \begin{aligned}
\Psi''\left(x\right) &= \frac{ \frac{\Psi\left(x+h\right) - \Psi\left(x\right)}{h} - \frac{\Psi\left(x\right) - \Psi\left(x-h\right)}{h}}{h}, \\ \\
&= \frac{\Psi\left(x+h\right) - 2\Psi\left(x\right) + \Psi\left(x-h\right)}{h^2}.
\end{aligned} $$

Thus the Schrodinger equation can be approximated by,

$$ \frac{-\hbar}{2m} \left( \frac{\Psi\left(x+h\right) - 2\Psi\left(x\right) + \Psi\left(x-h\right)}{h^2} 	\right) + V \left( x \right) \Psi \left( x \right) = E \Psi \left( x \right). $$

In matrix form, this discretized form of this equation at each step becomes,

$$
-\frac{\hbar^2}{2mh^2 }
\begin{pmatrix}
 -2 & 1 & 0 & ... & 0 \\
 1 & -2 & 1 & ... & 0 \\
\vdots &  &  & \ddots & \vdots \\
0 & & ... & -2 & 1 \\
0 & & ... & 1 & -2 \\
\end{pmatrix}
\begin{bmatrix}
\Psi \left( x_0 \right) \\
\Psi \left( x_1 \right) \\
\Psi \left( x_2 \right) \\
\vdots \\
\Psi \left( x_N \right) \\
\end{bmatrix}
+\begin{pmatrix}
V\left(x_0\right) & & & & \\
& V\left(x_1\right) & & & \\
& & V\left(x_2\right) & & \\
& & & \ddots & \\
& & & & V\left(x_N\right) \\
\end{pmatrix}
\begin{bmatrix}
\Psi \left( x_0 \right) \\
\Psi \left( x_1 \right) \\
\Psi \left( x_2 \right) \\
\vdots \\
\Psi \left( x_N \right) \\
\end{bmatrix} =
E
\begin{bmatrix}
\Psi \left( x_0 \right) \\
\Psi \left( x_1 \right) \\
\Psi \left( x_2 \right) \\
\vdots \\
\Psi \left( x_N \right) \\
\end{bmatrix},
$$

or,

$$
-\frac{\hbar^2}{2mh^2 } \left(
\begin{pmatrix}
 -2 & 1 & 0 & ... & 0 \\
 1 & -2 & 1 & ... & 0 \\
\vdots &  &  & \ddots & \vdots \\
0 & & ... & -2 & 1 \\
0 & & ... & 1 & -2 \\
\end{pmatrix}
+\begin{pmatrix}
 V\left(x_0\right) & & & & \\
 & V\left(x_1\right) & & & \\
 & & V\left(x_2\right) & & \\
 & & & \ddots & \\
 & & & & V\left(x_N\right) \\
\end{pmatrix}
\right)
\begin{bmatrix}
\Psi \left( x_0 \right) \\
\Psi \left( x_1 \right) \\
\Psi \left( x_2 \right) \\
\vdots \\
\Psi \left( x_N \right) \\
\end{bmatrix} =
E
\begin{bmatrix}
\Psi \left( x_0 \right) \\
\Psi \left( x_1 \right) \\
\Psi \left( x_2 \right) \\
\vdots \\
\Psi \left( x_N \right) \\
\end{bmatrix}.
$$

As discussed above, the left hand side of this equation is the Hamiltonian operator.  Since we are only interested in solutions that do not converge to zero at either end, the boundary conditions,

$$ \begin{aligned}
\Psi \left(x_0\right) &= 0 \\
\Psi \left(x_N\right) &= 0,
\end{aligned} $$

are applied.  The eigenstates of this system can therefore be found by diagonalizing the Hamiltonian operator, where the eigenvalues are the energies and the eigenvectors are the wavefunctions.
