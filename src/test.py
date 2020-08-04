#!/usr/bin/env python

import simulator as sim 
import geometry as geo 
from units import *
import ROOT as root
import read_data as rd
import physics
#https://indico.cern.ch/event/940711/
root.gStyle.SetPalette(55)
#EVENT_WRITE_DIR = "/global/project/projectdirs/dune/users/sgberg/sim/output/displays/"

zero = geo.point_3d(0., 0., 0.)

draw = False
events = [2]

def __main__():

	geometry = sim.default_geo()

	data_file = root.TFile.Open("~/hex/berk/nersc/data/jackBig.root")

	tree = data_file.Get("tree")
	
	dummy = rd.supervisor()
	d = rd.data(tree, dummy)
	d.stack()
	for k, event_number in enumerate(events):
		print("SIMULATING EVENT: " + str(event_number) )
		simul = sim.simulator(geometry)
		d = rd.data(tree, simul.supervisor)
		p, e, vtx = d.get_event_w_shift(event_number, module_box._xwalls, module_box._ywalls, module_box._zwalls )
		

		geom = simul.project_event(p, e, vtx)
		simul.get_stats()
		simul.supervisor.write()



__main__()
#







