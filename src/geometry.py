import numpy as np 
import ROOT as root 
from units import *
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cmx
from mpl_toolkits.mplot3d import Axes3D
EVENT_WRITE_DIR = "/home/stephen/hex/berk/sim/output/events/"

class point_3d:

	def __init__(self, x, y, z):
		self._x = x 
		self._y = y 
		self._z = z
		self._orig_x = 0.
		self._orig_y = 0.
		self._orig_z = 0.
		
	def shift_origion(self, x0, y0, z0): #default 0,0,0,
		self._orig_x = x0
		self._orig_y = y0
		self._orig_z = z0
	

	def get_xyz(self):
		return (self._x - self._orig_x, self._y- self._orig_y, self._z - self._orig_z)

	def get_rtp(self):
		x, y, z = (self._x - self._orig_x, self._y- self._orig_y, self._z - self._orig_z)
		r = np.sqrt(x*x + y*y + z*z)
		if r <= 0.0:
			theta = 0
		else:
			theta = np.arccos(z / r)

		if x == 0:
			phi = pi/2.0
		else:
			phi = np.arctan(y/x)

		return (r, theta, phi)

	def print(self):
		print("(" + str(self._x ) + ", " + str(self._y) + ", " + str(self._z) + ")" )

class wall:

	def __init__(self, th2, axis_position):

		self._phys = th2
		self._fixed_index = -1
		self._position = -1.0
		self._nx = th2.GetNbinsX()
		self._ny = th2.GetNbinsY()

		for k, val in enumerate(axis_position):
			if not val == 0:
				self._fixed_index = k
				self._position = val

	def will_hit(self, x, y, z, vector):
		
		if vector[self._fixed_index] == 0:
			return -1, []

		pt = [x, y, z]
		perp_distance = self._position - pt[self._fixed_index]

	
		if perp_distance/vector[self._fixed_index] < 0.0:
			return -1, []
		_bin = -1
		_bin_x = 0
		_bin_y = 0

		if self._fixed_index == 0:
			a = vector[0]
			xn = self._position
			yn = y + vector[1]*perp_distance/a 
			zn = z + vector[2]*perp_distance/a
			_bin_x = self._phys.GetXaxis().FindFixBin(yn)
			_bin_y = self._phys.GetYaxis().FindFixBin(zn)

		if self._fixed_index == 1:
			a = vector[1]
			xn = x + vector[0]*perp_distance/a 
			zn = z + vector[2]*perp_distance/a
			yn = self._position
			_bin_x = self._phys.GetXaxis().FindFixBin(xn)
			_bin_y = self._phys.GetYaxis().FindFixBin(zn)

		if self._fixed_index == 2:
			a = vector[2]
			yn = y + vector[1]*perp_distance/a 
			xn = x + vector[0]*perp_distance/a
			zn = self._position
			_bin_x = self._phys.GetXaxis().FindFixBin(xn)
			_bin_y = self._phys.GetYaxis().FindFixBin(yn)

		if _bin_x <= 0 or _bin_y <= 0:
			return -1, []
		if _bin_x > self._nx or _bin_y > self._ny:
			return -1, []

		_bin = self._phys.GetBin(_bin_x, _bin_y)

		return _bin, [xn, yn, zn] 
			

	def n_points(self):
		return self._phys.GetNbinsX(), self._phys.GetNbinsY()

	def get_point(self, j, l):
		point = []
		p1 = self._phys.GetXaxis().GetBinCenter(j)
		p2 = self._phys.GetYaxis().GetBinCenter(l)

		if (self._fixed_index == 0):
			point = [self._position, p1, p2]
		if (self._fixed_index == 1):
			point = [p1, self._position, p2]
		if (self._fixed_index == 2):
			point = [p1, p2, self._position]
		
		return (point[0], point[1], point[2])
		

class box:

	_name = ""
	_walls = []
	_NBINSX = int(1000/3)
	_NBINSY = int(3540/3)
	_NBINSZ = int(1000/3)

	def __init__(self, dimx, dimy, dimz, origion=point_3d(0.0, 0.0, 0.0), name="box" ):
		self._name = name
		self._dimx = dimx
		self._dimy = dimy 
		self._dimz = dimz
		self._origion = origion
		xw, yw, zw = self._origion.get_xyz()
		self._xwalls = []
		self._ywalls = []
		self._zwalls = []
		self.compute_walls()


	def compute_walls(self):
		xw, yw, zw = self._origion.get_xyz()
		self._xwalls = [xw - self._dimx/2.0, xw + self._dimx/2.0]
		self._ywalls = [yw - self._dimy/2.0, yw + self._dimy/2.0]
		self._zwalls = [zw - self._dimz/2.0, zw + self._dimz/2.0]

		self._walls = []
		self._walls.append(wall(root.TH2D(self._name + "_wallx_low" , "x=" + str(self._xwalls[0]), self._NBINSY, self._ywalls[0], self._ywalls[1], self._NBINSZ, self._zwalls[0], self._zwalls[1]) , [self._xwalls[0], 0, 0]))
		self._walls.append(wall(root.TH2D(self._name + "_wallx_high", "x=" + str(self._xwalls[1]), self._NBINSY, self._ywalls[0], self._ywalls[1], self._NBINSZ, self._zwalls[0], self._zwalls[1]) , [self._xwalls[1], 0, 0]))
		self._walls.append(wall(root.TH2D(self._name + "_wally_low" , "y=" + str(self._ywalls[0]), self._NBINSX, self._xwalls[0], self._xwalls[1], self._NBINSZ, self._zwalls[0], self._zwalls[1]) , [0, self._ywalls[0], 0]))
		self._walls.append(wall(root.TH2D(self._name + "_wally_high", "y=" + str(self._ywalls[1]), self._NBINSX, self._xwalls[0], self._xwalls[1], self._NBINSZ, self._zwalls[0], self._zwalls[1]) , [0, self._ywalls[1], 0]))
		self._walls.append(wall(root.TH2D(self._name + "_wallz_low" , "z=" + str(self._zwalls[0]), self._NBINSX, self._xwalls[0], self._xwalls[1], self._NBINSY, self._ywalls[0], self._ywalls[1]) , [0, 0, self._zwalls[0]]))
		self._walls.append(wall(root.TH2D(self._name + "_wallz_high", "z=" + str(self._zwalls[1]), self._NBINSX, self._xwalls[0], self._xwalls[1], self._NBINSY, self._ywalls[0], self._ywalls[1]) , [0, 0, self._zwalls[1]]))


	def is_inside(self, xp, yp, zp):
		xw, yw, zw = self._xwalls, self._ywalls, self._zwalls
		return ( (xp > xw[0] and xp < xw[1]) and (yp > yw[0] and yp < yw[1]) and (zp > zw[0] and zp < zw[1]))

	def set_origion(self, x, y, z):
		self._origion = point_3d(x, y, z)
		self.compute_walls() 

	def boundary(self): #returns wall array, normal vector to each wall , normal vector of boundary (arrays with corresponding indices)
		
		x_normal =   [-1.0, 0.0, 0.0]
		neg_x_norm = [ 1.0, 0.0, 0.0]
		y_normal =   [0.0, -1.0, 0.0]
		neg_y_norm = [0.0,  1.0, 0.0]
		z_normal =   [0.0, 0.0, -1.0]
		neg_z_norm = [0.0, 0.0,  1.0]

		normals = [x_normal, neg_x_norm, y_normal, neg_y_norm, z_normal, neg_z_norm]

		return self._walls, normals

	def fill(self, point):
		#point.print()
		return
	

class geometry:

	_E_field = 0.0

	def __init__(self, shapes, origion):
		self._shapes = shapes 
		self._origion = origion

	def is_inside(self, x, y, z):
		for k in range(len(self._shapes)):
			if self._shapes[k].is_inside(x, y, z):
				return k
		return -1

	def set_E(self, E):
		_E_field = E




def scatter3d(x,y,z, cs, tag, title="vox x,y,z vs. e", sym_z=True, colorsMap='jet'):
    cm = plt.get_cmap(colorsMap)
    cNorm = matplotlib.colors.Normalize(vmin=min(cs), vmax=max(cs))
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.scatter(x, y, z, s=1, c=scalarMap.to_rgba(cs))

    ax.set_xlim(-357.0, 357.0)
    ax.set_ylim(-150.0, 150.0)
    if sym_z:
    	ax.set_zlim(-255.0, 255.0)
    else:
    	ax.set_zlim(0.0, 510.0)


    ax.set_xlabel('X[cm]')
    ax.set_ylabel('Y[cm]')
    ax.set_zlabel('Z[cm]')
    ax.set_title(title)
    scalarMap.set_array(cs)
    fig.colorbar(scalarMap)
    plt.savefig("../plots/" + str(tag) + ".png")
    plt.close()

#geometric distance
def distance(pt1, pt2):
	return np.sqrt( (pt1[0]-pt2[0])*(pt1[0]-pt2[0]) + (pt1[1]-pt2[1])*(pt1[1]-pt2[1]) + (pt1[2]-pt2[2])*(pt1[2]-pt2[2]) )

def unit_vector(theta, phi):
	return [np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta)]

def rotate(vector, angle, axis):
	x, y, z = vector[0], vector[1], vector[2]
	a, b, c = axis[0], axis[1], axis[2]
	_cos = np.cos(angle)
	_sin = np.sin(angle)
	
	v1 = x*(_cos + a*a*(1-_cos)) + y*(a*b*(1-_cos) - c*_sin) + z*(a*c*(1-_cos) + b*_sin)
	v2 = x*(a*b*(1-_cos) + c*_sin) + y*(_cos + b*b*(1-_cos)) + z*(b*c*(1-_cos) - a*_sin)
	v3 = x*(a*c*(1-_cos) - b*_sin) + y*(b*c*(1-_cos) + a*_sin) + z*(_cos + c*c*(1-_cos))

	return [v1, v2, v3]


def scatter3d_module(x,y,z, cs, tag, title="vox x,y,z vs. e", sym_z=True, colorsMap='jet'):
    cm = plt.get_cmap(colorsMap)
    cNorm = matplotlib.colors.Normalize(vmin=min(cs), vmax=max(cs))
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.scatter(x, y, z, s=1, c=scalarMap.to_rgba(cs))

    ax.set_xlim(-50.0, 50.0)
    ax.set_ylim(-150.0, 150.0)
    ax.set_zlim(-50.0, 50.0)

    ax.set_xlabel('X[cm]')
    ax.set_ylabel('Y[cm]')
    ax.set_zlabel('Z[cm]')
    ax.set_title(title)
    scalarMap.set_array(cs)
    fig.colorbar(scalarMap)
    plt.savefig(EVENT_WRITE_DIR + "modules/event_" + str(tag) + ".png")
    plt.close()