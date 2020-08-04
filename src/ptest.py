import ROOT as root 
import numpy as np
import light_learner as ll
import pileup as pu
import matplotlib.pyplot as plt

def import_data(n):
	#print("../output/rands/rand_event" + str(n) + ".root")
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
	ylow = file.Get(prefix + "y_low")
	yhigh = file.Get(prefix + "y_high")
	zlow = file.Get(prefix + "z_low")
	zhigh = file.Get(prefix + "z_high")

	pts,  deps = pu.track_points([vtx_x, vtx_y, vtx_z], pu.unit_vector(theta, phi) )
	vals = pts[int(len(pts)/2)]

	walls = [xlow, xhigh, ylow, yhigh, zlow, zhigh]

	total_int = 0.0
	for k, wall in enumerate(walls):
		nx = wall.GetNbinsX()
		ny = wall.GetNbinsY()
		wall.Rebin2D(int(nx/25), int(ny/25))
		nx = wall.GetNbinsX()
		ny = wall.GetNbinsY()
		total_int += wall.Integral()

	x_data = []

	for k, wall in enumerate(walls):
	#	if wall.Integral() > 1.0:
		#	wall.Scale(1.0/total_int)
		for i in range(wall.GetNbinsX()):
			for j in range(wall.GetNbinsY()):
				x_data.append(wall.GetBinContent(i, j))
	x_data.append(total_int)


	file.Close()
	return x_data, vals

def get_truth(n):
	#print("../output/rands/rand_event" + str(n) + ".root")
	file = root.TFile.Open("../output/rands/rand_event" + str(n) + ".root")
	truth = file.Get("_truth")
	vtx_x = truth.GetBinContent(1)
	vtx_y = truth.GetBinContent(2)
	vtx_z = truth.GetBinContent(3)
	theta = truth.GetBinContent(4)
	phi = truth.GetBinContent(5)

	file.Close()
	return [vtx_x, vtx_y, vtx_z, theta, phi]

def get_center(n):
	file = root.TFile.Open("../output/rands/rand_event" + str(n) + ".root")
	truth = file.Get("_truth")
	vtx_x = truth.GetBinContent(1)
	vtx_y = truth.GetBinContent(2)
	vtx_z = truth.GetBinContent(3)
	theta = truth.GetBinContent(4)
	phi = truth.GetBinContent(5)
	pts,  deps = pu.track_points([vtx_x, vtx_y, vtx_z], pu.unit_vector(theta, phi) )
	vals = pts[int(len(pts)/2)]
	file.Close()
	return vals



fx, fy = import_data(87)

#inputs = np.ndarray((100, len(fx)))
#targets = np.ndarray((100, len(fy)))
inputs = []
targets = []

for k in range(350):
	x, y = import_data(k)
	inputs.append(x)
	targets.append(y)
	#for j in range(len(x)):
	#	inputs[k][j] = x[j]
	#for j in range(len(y)):
	#	targets[k][j] = y[j]


learner = ll.light_learner()
learner.train(inputs, targets)

for j in range(12):
	x, y1 = import_data(j+350)
	y = get_truth(j+350)
	print("truth:")
	print(y)
	out = learner._regressor.predict([x])

	pts, deps = pu.track_points([y[0], y[1], y[2]], pu.unit_vector(y[3], y[4]))
	print("guess:")
	outs = out[0]
	print(out)

	x, y, z = pu.unzip(pts)

	
	pu.scatter3d_module(x,y,z, deps, 10000+j, "truth vs. reconstructed center point", outs)


plot = root.TH1D("dist", "distance from correct center", 100, 0., 100.0)

for k in range(360):
	x, y1 = import_data(k)
	c = get_center(k)
	out = learner._regressor.predict([x])
	outs = out[0]
	x, y, z = outs[0] - c[0], outs[1]- c[1], outs[2] - c[2]
	dist = np.sqrt(x*x + y*y + z*z)
	plot.Fill(dist)

canv = root.TCanvas("c1", "c1")
plot.GetXaxis().SetTitle("distance [cm]")
plot.GetXaxis().CenterTitle()
plot.Draw()
canv.Print("dist.pdf", ".pdf" )




