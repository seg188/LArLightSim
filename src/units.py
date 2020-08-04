import numpy as np 

##DEFINING UNITS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#EVERYWHERE IN THIS SIMULATION, REFERENCE UNITS ARE MEV, CM, SECOND!!!!


#ENERGY
MeV = 1.00
GeV = np.power(10.0, 3.0)
keV = np.power(10.0, -3.0)
eV = np.power(10.0, -6.0)
kV = 1.00

#LENGTH

cm = 1.00
mm = np.power(10.0, -3.0)

#time

ns = 1.00
us = np.power(10.0, 3.0)

cc = 29.9*cm/ns #speed of light

pi = 3.1415