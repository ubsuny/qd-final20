from pylab import *
from scipy.integrate import simps
from SchrodingerSolver_user import Finite_Well_Solver

N_test = 1000
Vo_test = 1.515 #[eV]
L_test = 5*10**(-9) #[m]

def test_normalization(N, Vo, L, which_eig):
	b = 2*L #[m]
	x = linspace(-b, b, N)      # x-axis
	en = linspace(0, Vo, 100)  # vector of energies where we look for the stable states
	m = 0.51099895000*10**6 # [eV/c^2]
	hbar = 6.582119569*10**(-16) # Planck constant [eV*s]

	s = Finite_Well_Solver(x, b, L, Vo, en, m, hbar)

	eigvals, eigfuns = s.finite_diff(which_eig+1)

	assert simps(eigfuns[:,which_eig]**2) - 1.0 < 0.1


test_normalization(N_test, Vo_test, L_test, 1)
