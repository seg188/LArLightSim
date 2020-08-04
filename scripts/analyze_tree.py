import numpy as np 
import ROOT as root 
from units import *
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cmx
from mpl_toolkits.mplot3d import Axes3D
import physics

def scatter3d(x,y,z, cs, tag, colorsMap='jet'):
    cm = plt.get_cmap(colorsMap)
    cNorm = matplotlib.colors.Normalize(vmin=min(cs), vmax=max(cs))
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.scatter(x, y, z, s=1, c=scalarMap.to_rgba(cs))

    ax.set_xlim(-357.0, 357.0)
    ax.set_ylim(-150.0, 150.0)
    ax.set_zlim(-255.0, 255.0)

    ax.set_xlabel('X[cm]')
    ax.set_ylabel('Y[cm]')
    ax.set_zlabel('Z[cm]')
    ax.set_title("vox_x_y_z vs. photon yield")
    scalarMap.set_array(cs)
    fig.colorbar(scalarMap)
    plt.savefig("plots/event_" + str(tag) + ".pdf")
    plt.close()


def analyze_event(tree, event_n):
	tree.GetEntry(event_n)
	inter = physics.interaction(tree.fsPDG, tree.fsKE)
	

def analyze(tr):
	for k in range(tr.GetEntries()):
		analyze_event(tr, k)
	

def _____main_____():
	file_name  = "~/hex/berk/nersc/data/stephenBig.root"
	file = root.TFile.Open(file_name)

	tree = file.Get("tree")

	analyze(tree)

_____main_____()