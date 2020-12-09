from pylab import *
from scipy.integrate import odeint
from scipy.optimize import brentq

class Finite_Well_Solver:

	def __init__(self, x, b, L, Vo, en, psi, psi0):
		self.x = x
		self.b = b
		self.L = L
		self.Vo = Vo
		self.en = en
		self.psi = psi
		self.psi0 = psi0


	def root_bisection(self, f, x1, x2, accuracy=1.0e-6, max_steps=1000, root_debug=False):
	    """Return root of f(x) in bracketed by x1, x2 with specified accuracy.
	    Assumes that f(x) changes sign once in the bracketed interval.
	    Uses bisection root-finding algorithm.
	    """
	    iterations = []
	    f1 = f(x1)
	    f2 = f(x2)
	    if f1 * f2 > 0.0:
	        raise Exception("BISECTION METHOD: f(x1) * f(x2) > 0.0")
	    x_mid = (x1 + x2) / 2.0
	    f_mid = f(x_mid)
	    dx = x2 - x1
	    step = 0
	    if root_debug:
	        iterations = []
	        root_print_header("Bisection Search", accuracy)
	        root_print_step(step, x_mid, dx, f_mid)
	        iterations.append([x_mid,f_mid])
	    while abs(dx) > accuracy:
	        if f_mid == 0.0:
	            dx = 0.0
	        else:
	            if f1 * f_mid > 0:
	                x1 = x_mid
	                f1 = f_mid
	            else:
	                x2 = x_mid
	                f2 = f_mid
	            x_mid = (x1 + x2) / 2.0
	            f_mid = f(x_mid)
	            dx = x2 - x1
	        step += 1
	        if step > max_steps:
	            warning = "Too many steps (" + repr(step) + ") in root_bisection"
	            raise Exception(warning)
	        if root_debug:
	            root_print_step(step, x_mid, dx, f_mid)
	            iterations.append([x_mid,f_mid])
	    return x_mid,np.array(iterations)

	def V_func(self, inp):
		"""
    	Potential function in the finite square well. Width is L and value is global variable Vo
		"""
		if abs(inp) > self.L:
			return self.Vo
		else:
			return 0

	def SE(self, psi_inp, x_inp):
		"""
		Returns derivatives for the 1D schrodinger equation.
		Requires global value E to be set somewhere.  State0 is
		first derivative of the wave function psi, and state1 is
		its second derivative"""
		state0 = psi_inp[1]
		state1 = 2.0*(self.V_func(x_inp) - self.E)*psi_inp[0]
		return array([state0, state1])


	def Wave_function(self, energy):
		"""
		Calculates wave function psi for the given value of
		energy E and returns value at point b
		"""
		self.E = energy
		self.psi = odeint(self.SE, self.psi0, self.x)
		return self.psi[-1,0]

	def Wave_function_full(self, energy):
		"""
		Calculates wave function psi for the given value of
		energy E"""
		self.E = energy
		self.psi = odeint(self.SE, self.psi0, self.x)
		return self.psi


	def find_all_zeroes(self, x, y):
		"""
		Gives all zeroes in y = f(x)
		"""
		all_zeroes = []
		s = sign(y)
		for i in range(len(y)-1):
			if s[i]+s[i+1] == 0:
				zero = self.root_bisection(self.Wave_function, x[i], x[i+1])
				all_zeroes.append(zero)
		return all_zeroes

	def find_analytic_energies(self):
		"""
		Calculates Energy values for the finite square well using analytical
		model (Griffiths, Introduction to Quantum Mechanics, page 62.)
		"""
		z = sqrt(2*self.en)
		z0 = sqrt(2*self.Vo)
		z_zeroes = []
		z_ret = []
		f_sym = lambda z: tan(z)-sqrt((z0/z)**2-1)      # Formula 2.138, symmetrical case
		f_asym = lambda z: -1/tan(z)-sqrt((z0/z)**2-1)  # Formula 2.138, antisymmetrical case

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
			print("%.4f" %(z_zeroes[i]**2/2))
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
			print("%0.4f" % (z_zeroes_a[i]**2/2))


		return array([z_ret, z_ret_a])
