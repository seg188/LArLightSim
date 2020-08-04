import numpy as np 
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cmx
from mpl_toolkits.mplot3d import Axes3D
from units import *
import simulator as sim 
import geometry as geo 
import ROOT as root
import read_data as rd
import physics

EVENT_WRITE_DIR = "../plots/rands/"
OUT_DIR = "../output/rands/"
x_range = [-50.0, 50.0]#cm
y_range = [-150.0, 150.0]#cm
z_range = [-50.0, 50.0]#cm

POINT_SPACING = 5*mm


def scatter3d_module(x,y,z, cs, tag, title="track", vtx=[], sym_z=True, colorsMap='jet'):
    cm = plt.get_cmap(colorsMap)
    cNorm = matplotlib.colors.Normalize(vmin=min(cs), vmax=max(cs))
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.scatter(x, y, z, s=1, c=scalarMap.to_rgba(cs))

    ax.set_xlim(x_range[0], x_range[1])
    ax.set_ylim(y_range[0], y_range[1])
    ax.set_zlim(z_range[0], z_range[1])

    ax.set_xlabel('X[cm]')
    ax.set_ylabel('Y[cm]')
    ax.set_zlabel('Z[cm]')
    ax.set_title(title)
    scalarMap.set_array(cs)
    fig.colorbar(scalarMap)
    if len(vtx) == 3:
    	ax.scatter([vtx[0]], [vtx[1]], [vtx[2]], marker="X")

    plt.savefig(EVENT_WRITE_DIR + "rand_event_" + str(tag) + ".png")
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

def random_vertex(): #choose random vertex in middle 50% in all directions
	rands = np.random.rand(3)
	x = x_range[0] + (rands[0])*(x_range[1] - x_range[0])
	y = y_range[0] + (rands[1])*(y_range[1] - y_range[0])
	z = z_range[0] + (rands[2])*(z_range[1] - z_range[0])

	return [x, y, z]

def random_direction(): #choose random direction for track
	rands = np.random.rand(2)
	theta = rands[0]*pi
	phi   = rands[1]*pi
	return theta, phi

def inside_x(x):
	return (x > x_range[0] and x < x_range[1])
def inside_y(y):
	return (y > y_range[0] and y < y_range[1])
def inside_z(z):
	return (z > z_range[0] and z < z_range[1])

def inside_module(x, y, z):
	return (inside_x(x) and inside_y(y) and inside_z(z))

def track_points(vertex, direction):
	x0, y0, z0 = vertex[0], vertex[1], vertex[2]
	pts = [  ]
	deps = [  ]
	while inside_module(x0, y0, z0):
		pts.append([x0, y0, z0])
		deps.append(100)

		x0 = x0 + POINT_SPACING*direction[0]
		y0 = y0 + POINT_SPACING*direction[1]
		z0 = z0 + POINT_SPACING*direction[2]
		
	return pts, deps

def unzip(points):
	x = []
	y = []
	z = []
	for k, pt in enumerate(points):
		x.append(pt[0])
		y.append(pt[1])
		z.append(pt[2])
	return x, y, z

def make_random(k):
	vtx = random_vertex()
	theta, phi = random_direction()
	direct = unit_vector(theta, phi)
	pts, deps = track_points(vtx, direct)
	print(len(deps))
	x, y, z = unzip(pts)
	scatter3d_module(x, y, z, deps, k)

	zero = geo.point_3d(0., 0., 0.)
	#box = geo.box(714*cm, 300*cm, 510*cm, zero, "ar_cube")
	module_box = geo.box(100.0*cm, 300.0*cm, 100*cm, zero, "rand" + str(k) + "_box")
	geometry = geo.geometry([module_box], zero)
	simul = sim.simulator(geometry)
	simul.supervisor._outfile = OUT_DIR + "rand_event" + str(k) + ".root"

	geom = simul.project_event(pts, deps, vtx)
	simul.set_minimal()
	simul.supervisor._truth.SetBinContent(1, vtx[0])
	simul.supervisor._truth.SetBinContent(2, vtx[1])
	simul.supervisor._truth.SetBinContent(3, vtx[2])
	simul.supervisor._truth.SetBinContent(4, theta)
	simul.supervisor._truth.SetBinContent(5, phi)
	simul.supervisor.write()

def main():
	for k in range(1):
		make_random(281)
	

#main()