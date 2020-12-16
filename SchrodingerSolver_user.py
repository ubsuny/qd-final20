from pylab import *
from scipy.integrate import odeint
from scipy.optimize import brentq
import numpy as np

class Finite_Well_Solver:

	def __init__(self, x, b, L, Vo, en, m, hbar):
		self.x = x
		self.b = b
		self.L = L
		self.Vo = Vo
		self.en = en
		self.m = m
		self.hbar = hbar

	def V_func(self, inp):
		"""
		Potential function in the finite square well. Width is L and value is global variable Vo
		"""
		if abs(inp) > self.L:
			return self.Vo
		else:
			return 0


	def build_Hamiltonian(self):
		N = len(self.x)
		h = abs(self.x[1]-self.x[0])
		T = np.zeros((N-2)**2).reshape(N-2,N-2)
		for i in range(N-2):
			for j in range(N-2):
				if i==j:
					T[i,j] = -2
				elif np.abs(i-j)==1:
					T[i,j] = 1
				else:
					T[i,j] = 0
		#T = -T/(2*(h**2))
		T = -((self.hbar**2)/(2*self.m*(h**2))) * T
		V = np.zeros((N-2)**2).reshape(N-2,N-2)
		for i in range(N-2):
			for j in range(N-2):
				if i==j:
					V[i,j]= self.V_func(self.x[i+1])
				else:
					V[i,j]=0

		return T + V

	def finite_diff(self, num_eigs):
		H = self.build_Hamiltonian()

		val,vec=np.linalg.eig(H)
		eigval_ind = np.argsort(val)
		eigval_ind = eigval_ind[0:num_eigs]
		energies = val[eigval_ind]

		wavefunctions = np.zeros((len(val), num_eigs))
		for i in range(len(eigval_ind)):
			wavefunctions[:,i] = vec[:,eigval_ind[i]]

		zero_vec = np.zeros(num_eigs)
		wavefunctions = np.vstack((zero_vec,wavefunctions,zero_vec))
		wavefunctions[:,0] = abs(wavefunctions[:,0])

		return energies, wavefunctions

	def find_analytic_energies(self, sigfigs):
		"""
		Calculates Energy values for the finite square well using analytical
		model (Griffiths, Introduction to Quantum Mechanics, page 62.)
		"""
		z = (sqrt(2*self.m*(self.en))*self.L)/self.hbar
		z0 = (sqrt(2*self.m*self.Vo)*self.L)/self.hbar
		# z = (sqrt(2*self.en)*self.L)
		# z0 = (sqrt(2*self.Vo)*self.L)
		z_zeroes = []
		z_ret = []
		f_sym = lambda z: tan(z)-sqrt((z0/z)**2-1)      # Formula 2.156, symmetrical case
		f_asym = lambda z: -1/tan(z)-sqrt((z0/z)**2-1)  # Formula 2.156, antisymmetrical case

		# first find the zeroes for the symmetrical case
		s = sign(f_sym(z))
		for i in range(len(s)-1):   # find zeroes of this crazy function
			if s[i]+s[i+1] == 0:
				zero = brentq(f_sym, z[i], z[i+1])
				z_zeroes.append(zero)
		print("Energies from the analytical model are: ")
		print("(Symmetrical case)")
		for i in range(0, len(z_zeroes),2):   # discard z=(2n-1)pi/2 solutions cause that's where tan(z) is discontinuous
			z_ret.append(z_zeroes[i]**2/2)
			print("{:.{}f}".format(z_zeroes[i]**2/2, sigfigs))
		# Now for the asymmetrical
		z_zeroes_a = []
		z_ret_a = []
		s = sign(f_asym(z))
		for i in range(len(s)-1):   # find zeroes of this crazy function
			if s[i]+s[i+1] == 0:
				zero = brentq(f_asym, z[i], z[i+1])
				z_zeroes_a.append(zero)
		print("(Antisymmetrical case)")
		for i in range(0, len(z_zeroes_a),2):   # discard z=npi solutions cause that's where cot(z) is discontinuous
			z_ret_a.append(z_zeroes_a[i]**2/2)
			print("{:.{}f}".format(z_zeroes_a[i]**2/2, sigfigs))


		return array([z_ret, z_ret_a])
