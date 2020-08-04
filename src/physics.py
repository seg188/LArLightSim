from units import *
import bethebloch as bb
import numpy as np

gamma_e = 9.686*eV

def doke_birks(_E_d, dedx=2.19*MeV/cm):
	E_d = _E_d/10.0
	dokeBirks = [0.0, 0.0, 0.0]
	dokeBirks[0] = 0.07*np.power(E_d,-0.85)
	dokeBirks[2] = 0.00
	dokeBirks[1] = dokeBirks[0]/(1-dokeBirks[2]) #B=A/(1-C) (see paper)
	recombProb = (dokeBirks[0]*dedx)/(1.00+dokeBirks[1]*dedx) + dokeBirks[2]
	return recombProb
    

def light_yield(energy_dep, db_factor):
	exc_ratio = 0.25
	scint_yield = 1.0/(19.5*eV)
	num_quanta = energy_dep*scint_yield
	num_exciton = exc_ratio*num_quanta/(1.0 + exc_ratio)
	num_ions = num_quanta - num_exciton
	num_photons = num_exciton + db_factor*num_ions

	return num_photons

def rayleigh_differential_xs(theta):
	return 1.00 + np.cos(theta)*np.cos(theta)

def emission_spectrum(t):
	slow_ampl = 0.75
	fast_ampl = 0.25
	tau_fast = 5.0*ns
	tau_slow = 1.0*us

	return fast_ampl/tau_fast*np.exp(-1.0*t/tau_fast) + slow_ampl/tau_slow*np.exp(-1.0*t/tau_slow)

#######################################################################################################3
class interaction:

	def name(self, pdg):
		if pdg > 4000 or pdg < -4000:
			print(pdg)
		if pdg == 11:
			return "e-"
		if pdg == -11:
			return "e+"
		if pdg == 12:
			return "ve"
		if pdg == -12:
			return "_ve"
		if pdg == 13:
			return "mu-"
		if pdg == -13:
			return "mu+"
		if pdg == 14:
			return "vu"
		if pdg == -14:
			return "_vu"
		if pdg == 22:
			return "gamma"
		if pdg == 111:
			return "pi0"
		if pdg == 211:
			return "pi+"
		if pdg == -211:
			return "pi-"
		if pdg == 2112:
			return "n"
		if pdg == -2112:
			return "_n"
		if pdg == 2212:
			return "p"
		if pdg == -2212:
			return "_p"
		if pdg == 3112:
			return "sigma-"
		if pdg == 3222:
			return "sigma+"
		if pdg == 3122:
			return "lambda"
		if pdg == -3122:
			return "_lambda"
		if pdg == 311:
			return "K0"
		if pdg == -311:
			return "_K0"
		if pdg == 321:
			return "K+"
		if pdg == -321:
			return "K-"
		if pdg == 130:
			return "KL"
		return "unknown"

	_current = ""
	_class   = ""
	_flavor  = ""

	def __init__(self, pdgs, energies, voxes=[]):
		self.pdgs = pdgs
		self._es   = energies
		self._particles = []
		self._total_dep = 0.0
		for k, e in enumerate(voxes):
			self._total_dep = self._total_dep + e
		for k, pdg in enumerate(pdgs):
			self._particles.append(self.name(pdg))

		self._energy = 0.000
		for k, energy in enumerate(energies):
			self._energy = self._energy + energy

		self.classify()

	def has(self, name):
		for k, particle in enumerate(self._particles):
			if name == particle:
				return True
		return False

	def number_of(self, name):
		n = 0
		for k, particle in enumerate(self._particles):
			if name == particle:
				n = n + 1
		return n

	def n_hadrons(self):
		n = 0
		for k, particle in enumerate(self._particles):
			if particle == "neutron":
				n = n + 1
				continue
			if particle == "n":
				n = n + 1
				continue
			if particle == "_n":
				n = n + 1
				continue
			if particle == "p":
				n = n + 1
				continue
			if particle == "_p":
				n = n + 1
				continue
			if particle == "pi0":
				n = n + 1
				continue
			if particle == "pi-":
				n = n + 1
				continue
			if particle == "pi+":
				n = n + 1
				continue
			if particle == "sigma-":
				n = n + 1
				continue
			if particle == "sigma+":
				n = n + 1
				continue
			if particle == "lambda":
				n = n + 1
				continue
			if particle == "_lambda":
				n = n + 1
			if particle == "K0":
				n = n + 1
			if particle == "K+":
				n = n + 1
			if particle == "_K0":
				n = n + 1
			if particle == "KL":
				n = n + 1	
		return n

	def has_neutrino(self):
		return (self.has("vu") or self.has("_vu") or self.has("ve") or self.has("_ve"))
	def has_lepton(self):
		return (self.has("mu+") or self.has("mu-") or self.has("e+") or self.has("e-"))

	def n_non_nucleon_baryons(self):
		n = 0
		for k, particle in enumerate(self._particles):
			if particle == "sigma-":
				n = n + 1
				continue
			if particle == "sigma+":
				n = n + 1
				continue
			if particle == "lambda":
				n = n + 1
				continue
			if particle == "_lambda":
				n = n + 1
			if particle == "K0":
				n = n + 1
			if particle == "K+":
				n = n + 1
			if particle == "_K0":
				n = n + 1
			if particle == "KL":
				n = n + 1	
		return n

	def particle_string(self):
		string = "["
		for k, particle in enumerate(self._particles):
			string = string + " " + particle 
		string = string + " ]"
		return string

	def classify(self):
		current = "unkownn"
		_class  = "unknown"
		_flavor = "unknown"

		#GETTING FLAVOR OF INCIDENT MUON FOR NON-ELECTRON SCATTERING EVENTS

		l_u = self.number_of("mu-") + self.number_of("vu") - (self.number_of("mu+") + self.number_of("_vu"))
		l_e = self.number_of("e-") + self.number_of("ve") - (self.number_of("e+") + self.number_of("_ve"))



		if l_u == 1 and l_e == 0:
			_flavor = "MUON"
		if l_u == -1 and l_e == 0:
			_flavor = "ANTIMUON"
		if l_u == 0 and l_e == 1:
			_flavor = "ELECTRON"
		if l_u == 0 and l_e == -1:
			_flavor = "ANTIELECTRON"
		#######################################################################################################

		if _flavor == "MUON" and self.has("mu-"):
			current = "CHARGED_CURRENT"
		if _flavor == "ANTIMUON" and self.has("mu+"):
			current = "CHARGED_CURRENT"
		if _flavor == "ELECTRON" and self.has("e-"):
			current = "CHARGED_CURRENT"
		if _flavor == "ANTIELECTRON" and self.has("e+"):
			current = "CHARGED_CURRENT"

		if _flavor == "MUON" and self.has("vu"):
			current = "NEUTRAL_CURRENT"
		if _flavor == "ANTIMUON" and self.has("_vu"):
			current = "NEUTRAL_CURRENT"
		if _flavor == "ELECTRON" and self.has("ve"):
			current = "NEUTRAL_CURRENT"
		if _flavor == "ANTIELECTRON" and self.has("_ve"):
			current = "NEUTRAL_CURRENT"

		n_had = self.n_hadrons()
		n_res_bars = self.n_non_nucleon_baryons()
	

		if n_had == 1:
			if self.has("p") or self.has("n"):
				if self.has_lepton() and len(self._particles) == 2:
					_class = "QUASIELASTIC"
				if self.has_neutrino() and len(self._particles) == 2:
					_class = "QUASIELASTIC"

				if self.has("gamma") and n_res_bars == 0:
					_class = "INELASTIC_NUCLEON"

			if (self.has("pi0") or (self.has("pi+") or self.has("pi-"))) and len(self._particles) == 2:
				_class = "PION_PRODUCTION"

		if n_had >= 2:
			_class = "DEEP_INELASTIC"


		if n_had == 0:
			if self.has("mu-") and self.has("ve"):
				_class = "NEUTRINO_ELECTRON_SCATTERING"
				_flavor = "MUON"
				current = "CHARGED_CURRENT"
			if self.has("e-") and self.has("vu"):
				_class = "NEUTRINO_ELECTRON_SCATTERING"
				_flavor = "MUON"
				current = "NEUTRAL_CURRENT"


		if _class == "unknown":
			print(self._particles)

		self._class = _class
		self._flavor = _flavor
		self._current = current










######################################################################################################################################





   

    