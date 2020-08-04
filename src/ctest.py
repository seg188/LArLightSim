import ROOT as root 
import numpy as np
import light_learner as ll
import pileup as pu
import matplotlib.pyplot as plt

def import_data(n):
	print("../output/rands/rand_event" + str(n) + ".root")
	file = root.TFile.Open("../output/rands/rand_event" + str(n) + ".root")
	truth = file.Get("_truth")
	vtx_x = truth.GetBinContent(1)
	vtx_y = truth.GetBinContent(2)
	vtx_z = truth.GetBinContent(3)
	theta = truth.GetBinContent(4)
	phi = truth.GetBinContent(5)

	prefix = "rand" + str(n) + "_box_wall"
	xlow = file.Get(prefix + "x_low")
	xhigh = file.Get(prefix + "x_high")
	zlow = file.Get(prefix + "z_low")
	zhigh = file.Get(prefix + "z_high")

	#pts,  deps = pu.track_points([vtx_x, vtx_y, vtx_z], pu.unit_vector(theta, phi) )
	#vals = pts[int(len(pts)/2)]
	vals = [vtx_x, vtx_y, vtx_z]
	walls = [xlow, xhigh, zlow, zhigh]

	x_data = []

	total_int = 0.0
	for k, wall in enumerate(walls):
		wall.Rebin2D(5, 5)
		total_int += wall.Integral()


	for k, wall in enumerate(walls):
		if wall.Integral() > 1.0:
			wall.Scale(1.0/total_int)

		for i in range(wall.GetNbinsX()):
			for j in range(wall.GetNbinsY()):
				x_data.append(wall.GetBinContent(i, j))


	file.Close()
	return x_data, vals

def get_truth(n):
	print("../output/rands/rand_event" + str(n) + ".root")
	file = root.TFile.Open("../output/rands/rand_event" + str(n) + ".root")
	truth = file.Get("_truth")
	vtx_x = truth.GetBinContent(1)
	vtx_y = truth.GetBinContent(2)
	vtx_z = truth.GetBinContent(3)
	theta = truth.GetBinContent(4)
	phi = truth.GetBinContent(5)

	file.Close()
	return [vtx_x, vtx_y, vtx_z, theta, phi]


fx, fy = import_data(87)

#inputs = np.ndarray((100, len(fx)))
#targets = np.ndarray((100, len(fy)))
inputs = []
targets = []

for k in range(300):
	x, y = import_data(k)
	inputs.append(x)
	targets.append(1.0)
	#for j in range(len(x)):
	#	inputs[k][j] = x[j]
	#for j in range(len(y)):
	#	targets[k][j] = y[j]



learner = ll.light_learner()
learner.train(inputs, targets)

for j in range(12):
	x, y1 = import_data(j+350)
	y = get_truth(j+350)
	print(y)
	out = learner._regressor.predict([x])

	pts, deps = pu.track_points([y[0], y[1], y[2]], pu.unit_vector(y[3], y[4]))

	outs = out[0]
	print(outs)

	x, y, z = pu.unzip(pts)

	
	pu.scatter3d_module(x,y,z, deps, 10000+j, "truth vs. reconstructed center point", outs)

