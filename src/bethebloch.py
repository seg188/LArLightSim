import numpy as np
import ROOT as root

eV = np.power(10.0, -6.0)
N_A = 6.022*np.power(10.0, 23.0)

r_e = 2.817*np.power(10.0, -13.0)
pi=3.1415
K = 0.307075 #MeV/mol cm^2
#MUON CHARGE NMBER
z = 1.0
#lAr Atomic number
Z = 18.0
A = 39.95
rho = 1.39

hbar_omega = 28.816*np.sqrt(rho*Z/A)*np.power(10.0, -6.0)

#density corecion to ionization loss
delta_bg = np.sqrt(rho*(Z/A))

mu_mass = 105.7 #MEV
m_e = 0.511

#mean excitation emergy
I = 215.0*np.power(10.0, -6.0) #NOTE!! IN EV

#maximum energy transfer to electron!!!

def gamma(E): #ONLY FOR MUON
	return E/mu_mass

def beta_2(E): #VALID ONLY FOR A MUON
	g = gamma(E)
	return 1.000 - 1.000/(g*g)

def W_m(E):
	val = 2*m_e*beta_2(E)*gamma(E)*gamma(E)
	val = val/(1.0 + 2.0*gamma(E)*m_e/mu_mass + np.power(m_e/mu_mass, 2.0))
	return -1.0*val

def delta(E):
	beta = np.sqrt(beta_2(E))
	gam = gamma(E)

	return np.log(beta*gam) - 0.50 + np.log(hbar_omega/I)

def p(E):
	return np.sqrt(E*E - mu_mass*mu_mass)
		
def dEdx(E):
	T_CUT = 2.0
	#val = (0.5*np.log(2*m_e*beta_2(E)*gamma(E)*gamma(E)*T_CUT/(I)) - beta_2(E)*(1.0+T_CUT/W_m(E))/2.0 + delta(E) )
	#val = rho*A*val*-1.0*K*z*z*(Z/A)/(beta_2(E))
	#factor = 10.0*rho*4.0*pi*N_A*r_e*r_e*m_e*z*z*(Z/A)*(1.0/beta_2(E))
	factor = 0.189/beta_2(E)
	val = factor*(np.log(2.0*m_e*gamma(E)*gamma(E)*beta_2(E)/I) - beta_2(E))# - delta(E))

	return val

def dEdx_p(p):
	E = np.sqrt(p*p + mu_mass*mu_mass)
	return dEdx(E)



