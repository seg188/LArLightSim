import numpy as np 
from geometry import *
import geometry as geo
from units import *
import physics
from physics import rayleigh_differential_xs
import sgrandom as rdm
from read_data import supervisor

##################################################################################
#DEFAULT GEOMETRY FOR SIMULATIONS
#JUST INSTANTIATE GEOMETRY AS default_geo() 
##################################################################################
def default_geo():
	module_box = geo.box(100.0*cm, 300.0*cm, 100*cm)
	geometry = geo.geometry([module_box], zero)
	return geometry

class simulator: #basic class to handle projections, etc
	_project_area = 1
	_t = 0.0*ns
	_use_rayleigh = True
	_scattering_length = 95*cm
	_std_max = 511.0 #ensures cuttoff < .1%
	_diff_xs = rdm.Distribution(rayleigh_differential_xs, 0.0, pi)
	_scatter_radius = rdm.Distribution(rdm.poisson(_scattering_length), 0.0, _std_max)
	_emission_spectrum = rdm.Distribution(physics.emission_spectrum, 0.0, 3000.0*ns) #high precission sampling of emission time spectrum

	_scattered_points = []
	_scattered_e      = []

	_nscatters = 0
	_total_distance = 0.0

	_time_offset = 0.0
	_efficiency = 0.200
		

	def __init__(self, geometry):	
		self._geometry = geometry
		self.supervisor = supervisor()


	def rayleigh_project(self, point, strength):
		print("rayleigh scattering: " + str(strength))

		#simulates isotropic photons individually, uses raleigh scattering model with _scattering_length defined above

	def uniform_sim(self, x0, y0, z0, scatter=False, unit_vect=[]):

		container_index = self._geometry.is_inside(x0, y0, z0)
		
		if container_index == -1:
			return -1

		#container that has point inside
		container = self._geometry._shapes[container_index]
		bwalls, bnormals = container.boundary()

		theta = 0.0
		phi = 0.0

		if (not scatter):

			theta = (np.random.rand(1)*pi)[0]
			phi  = (np.random.rand(1)*2*pi)[0]
			self.supervisor._photon_theta_phi.Fill(theta, phi)
			unit_v = unit_vector(theta, phi)

		else:
			unit_v = unit_vect

		#unit vector in direction of propagation
		
		hit_wall = False
		for k, t_wall in enumerate(bwalls):
			wall = t_wall._phys
			wh, pt = t_wall.will_hit(x0, y0, z0, unit_v)
			if wh > 0:

				hit_wall = True
			

				dist = distance(pt, [x0, y0, z0])
				scatter_r = (self._scatter_radius.sample(1))[0]

				if scatter_r < dist:
					#do scattering here
					self._nscatters = self._nscatters + 1 
					self._total_distance = self._total_distance + scatter_r
					##################################################
					self.supervisor._scatter_distance.Fill(scatter_r)
					#################################################
					xn, yn, zn = x0 + unit_v[0]*dist, y0 + unit_v[1]*dist, z0 + unit_v[2]*dist
					new_theta = theta + (self._diff_xs.sample(1))[0]
					to_rotate = unit_vector(new_theta, phi)

					new_phi   = (np.random.rand(1))[0]*2.0*pi
					new_vector = rotate(to_rotate, new_phi, unit_v)
					self.supervisor._scatter_theta.Fill( new_vector[0]*unit_v[0] + new_vector[1]*unit_v[1] + new_vector[2]*unit_v[2] )
					val = self.uniform_sim(xn, yn, zn, True, new_vector)
					return val

				else:
					self._total_distance = self._total_distance + dist
					self.supervisor._detection_time.Fill(self._total_distance / cc )
					self._total_distance = 0.0
					self.supervisor._n_scatters.Fill(self._nscatters)
					self._nscatters = 0
					wall.SetBinContent(wh, wall.GetBinContent(wh) + 1.0) 
					return 0

		if not hit_wall:
			print("no hits!!")

		return 0
		
		

#fast project uses solid angle approximation, mre of statistical approach
#results are unreliable, and it isnt actually faster than the uniform_sim function

	def fast_project(self, point, strength=1.0, is_rayleigh=False): 
		x0, y0, z0 = point.get_xyz()
		xc, yc, zc = point.get_xyz()
		container_index = self._geometry.is_inside(x0, y0, z0)
		_sum = 0.0
		if container_index == -1:
			return [], []

		container = self._geometry._shapes[container_index]

		bwalls, bnormals = container.boundary()

		coefs = []

		_scatterting_function = rdm.poisson(self._scattering_length)

		for k, t_wall in enumerate(bwalls):
			wall = t_wall._phys
			norm = bnormals[k]
			x0 = container._xwalls[0]
			nx, ny = t_wall.n_points()
			for i in range(nx):
				for j in range(ny):

					(x0, y0, z0) = t_wall.get_point(i, j)

					point.shift_origion(x0, y0, z0)
					r, theta, phi = point.get_rtp()
					
					#print(r)
					x1, y1, z1 = point.get_xyz()

	
					unit_r = [x1/r, y1/r, z1/r]
					

					coef = 0
					for l in range(len(unit_r)):
						coef = coef + unit_r[l] * norm[l]


					total_e = np.absolute(coef)*strength/(4.0*pi*r*r) #energy per 1 square cm
				
					
					non_scattered = total_e*_scatterting_function(r)
					if non_scattered/(10.0*eV) > 1.0:
						print(non_scattered/(10.0*eV))
					
					if self._use_rayleigh:
						scattered = total_e - non_scattered
						if scattered > 10*eV:
							nphoton = int((scattered/(10.0*eV)))
							rs = np.random.rand(nphoton)*r
							phi = np.random.rand(nphoton)*2*pi
							theta = self._diff_xs.sample(nphoton)
							for ij, ir in enumerate(rs):
								r = [ir * d for d in unit_r]
								new_point = point_3d(xc + r[0], yc + r[1], zc + r[2])
								d = ir*np.cos(theta[ij])
								#print(d)
								bin1 = int(np.cos(phi[ij])*d/10.0)
								bin2 = int(np.sin(phi[ij])*d/10.0)
								print(bin1, bin2)
								_bin = wall.GetBin(i+bin1, j+bin2)
								wall.SetBinContent(_bin, wall.GetBinContent(_bin) + 1.00)

					if non_scattered/(10.0*eV) > 0.50:
						wall.SetBinContent(i, j, wall.GetBinContent(i, j) + round(non_scattered/(10.0*eV)) )
					else:
						rando = np.random.rand(1)*0.50
						if rando < non_scattered/(10.0*eV):
							wall.SetBinContent(i, j, wall.GetBinContent(i, j) + 1.0 )

	def set_minimal(self):
		walls, ns = (self._geometry._shapes[0]).boundary()
		for k, wall in enumerate(walls):
			t_wall = wall._phys
			self.supervisor.add_extras([t_wall])

	def get_stats(self):
		
		walls, ns = (self._geometry._shapes[0]).boundary()

		for k, wall in enumerate(walls):
			t_wall = wall._phys
			print(t_wall.GetName())
			for i in range(t_wall.GetNbinsX()):
				for j in range(t_wall.GetNbinsY()):
					self.supervisor._n_photon_per_pad.Fill(t_wall.GetBinContent(i, j))

			self.supervisor.add_extras([t_wall])

	def project_event(self, points, strenghts, vtx):
		l = float(len(points))
		n_inside = 0

		def distance_to_vtx(x, y, z):
			xp, yp, zp = x-vtx[0], y-vtx[1], z-vtx[2]
			return np.sqrt(xp*xp + yp*yp + zp*zp)
	
		for k, pt in enumerate(points):
			nextpt = False
			print(float(k)/l*100)
			if not nextpt:
				for n in range(int(strenghts[k]*self._efficiency)):
					self._total_distance = distance_to_vtx(pt[0], pt[1], pt[2])
					self._total_distance = cc * (self._emission_spectrum.sample()[0] + self._time_offset)
					status = self.uniform_sim(pt[0], pt[1], pt[2])#, strenghts[k])
					if status == -1:
						nextpt = True

		
		return self._geometry

	def clear_event(self):
		self.supervisor = supervisor()
		

	


		

	






	



