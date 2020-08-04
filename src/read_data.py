import ROOT as root
import numpy as np
import physics
from units import *
import geometry as geo

WRITE_DIR = "~/hex/berk/sim/output/root/modules/"


class supervisor:
	_outfile = "sim_output.root"

	def __init__(self, tag=""):
		self._extras = []
		self._scatter_distance = root.TH1D(tag + "_scatter_distance", "scattering distance",     100,    0.0,  500.00)
		self._n_scatters       = root.TH1D(tag + "_n_scatters",    "n scatters per photon",      20,   -0.50,   19.50)
		self._n_photon_per_pad = root.TH1D(tag + "_n_photon_per_pad", "photons per pad (3 mm)^2",2000, -0.50, 1999.50)
		self._sd_photon_per_pad= root.TH1D(tag + "_sd_photon_per_pad","sd of photons per pad"   ,100,   0.0,   100.0 )
		self._photon_yield     = root.TH1D(tag + "_photon_yield",  "photons per MeV deposit",    1000, 25000.0, 55000.0)
		self._scatter_theta    = root.TH1D(tag + "_scatter_theta", "cos(theta_scatter)"       ,  100,   -1.20,   1.20 )
		self._detection_time   = root.TH1D(tag + "_detection_time", "detection time",            1000, -5.0*ns,   3000*ns )
		self._photon_theta_phi = root.TH2D(tag + "_photon_theta_phi","gen photon theta vs phi",  100,   0.0,    pi,     100, 0.0, 2*pi)
		self._truth            = root.TH1D(tag + "_truth", "x,y,z,theta,phi",                     5,    -0.50, 4.50)

	def write(self):
		file = root.TFile(self._outfile, "RECREATE")
		file.cd()
		self._scatter_distance.Write()
		self._n_scatters.Write()
		self._n_photon_per_pad.Write()
		self._sd_photon_per_pad.Write()
		self._photon_yield.Write()
		self._scatter_theta.Write()
		self._photon_theta_phi.Write()
		self._detection_time.Write()
		self._truth.Write()
		for k, plot in enumerate(self._extras):
			plot.Write()
		file.Close()

	def add_extras(self, extra_plots):
		for k, plot in enumerate(extra_plots):
			self._extras.append(plot)

class data:

	_E_FIELD = 0.50*kV/cm

	def __init__(self, tree, supervisor):
		self._tree = tree 
		self.supervisor = supervisor

	def get_event(self, n):
		self._tree.GetEntry(n)
		self.supervisor._outfile = 	WRITE_DIR + "sim_out_event" + str(n) + ".root"
		deps = []
		points = []

		print(self._tree.nvox)
		for k in range(self._tree.nvox):
			db = physics.doke_birks(self._E_FIELD, self._tree.vox_dedx[k]/10.0)
			
			points.append([self._tree.voxx[k], self._tree.voxy[k], self._tree.voxz[k] - 255.0 ]) #shifting z to be symmetric

			edep = self._tree.voxe[k]*MeV

			_ly = physics.light_yield(edep, db)
			self.supervisor._photon_yield.Fill(_ly/edep)

			deps.append( _ly )

		vtx = [self._tree.vtxx, self._tree.vtxy, self._tree.vtxz] 

		return points, deps, vtx

	def get_event_w_shift(self, n, x_range, y_range, z_range):
		self._tree.GetEntry(n)
		self.supervisor._outfile = 	WRITE_DIR + "sim_out_event" + str(n) + ".root"
		deps = []
		points = []

		new_vtx_rands = np.random.rand(3)*0.50 
		vtx_x = x_range[0] + (new_vtx_rands[0] + 0.25)*(x_range[1]-x_range[0])
		vtx_y = y_range[0] + (new_vtx_rands[1] + 0.25)*(y_range[1]-y_range[0])
		vtx_z = z_range[0] + (new_vtx_rands[2] + 0.25)*(z_range[1]-z_range[0])

		shift = [vtx_x - self._tree.vtxx, vtx_y - self._tree.vtxy, vtx_z - self._tree.vtxz -5.0]

		print(self._tree.nvox)
		for k in range(self._tree.nvox):
			db = physics.doke_birks(self._E_FIELD, self._tree.vox_dedx[k]/10.0)
			
			points.append([self._tree.voxx[k] + shift[0], self._tree.voxy[k] + shift[1], self._tree.voxz[k] + shift[2] ]) #shifting z to be symmetric

			edep = self._tree.voxe[k]*MeV

			_ly = physics.light_yield(edep, db)
			self.supervisor._photon_yield.Fill(_ly/edep)

			deps.append( _ly )

		vtx = [vtx_x, vtx_y, vtx_z]

		return points, deps, vtx
		


	def study(self):
		plot = root.TGraph(1000)
		for k in range(1000):
			dedx = 0.50 + 4.0*float(k)/1000.0
			db = physics.doke_birks(self._E_FIELD, dedx)
			plot.SetPoint(plot.GetN(), dedx, physics.light_yield(1.00*MeV, db))

		plot.SetTitle("dE/dx Dependence of Light Yield at E=0.5 kV/cm")
		plot.GetXaxis().SetTitle("dE/dx [MeV/cm]")
		plot.GetXaxis().CenterTitle()

		plot.GetYaxis().SetTitle("Photons per MeV")
		plot.GetYaxis().CenterTitle()

		canv = root.TCanvas("canv", "canv")

		plot.SetLineColor(1)
		plot.SetLineWidth(2)

		plot.Draw("AC")
		canv.Print("ly_dedx_dep.pdf", ".pdf")
		return

	flavors = ["ELECTRON", "ANTIMUON", "ANTIELECTRON", "MUON"]
	currents = ["NEUTRAL_CURRENT", "CHARGED_CURRENT"]
	classes = [ "PION_PRODUCTION","INELASTIC_NUCLEON","NEUTRINO_ELECTRON_SCATTERING", "QUASIELASTIC", "DEEP_INELASTIC"]

	def draw_characteristic_events(self):
		for i, curr in enumerate(self.currents):
			for j, _class in enumerate(self.classes):
				low_e = False
				high_e = False

				for k in range(self._tree.GetEntries()):
					self._tree.GetEntry(k)
					inter = physics.interaction(self._tree.fsPDG, self._tree.fsKE)

					if ((inter._class == _class) and (inter._current == curr)) and self._tree.nvox > 50:
						if (not low_e) and inter._energy < 10*GeV:
							low_e = True
							geo.scatter3d(self._tree.voxx, self._tree.voxy, self._tree.voxz, self._tree.voxe, "event_" + str(k), inter.particle_string() + " KE=" + str(round(inter._energy/GeV) ), False)
						if (not high_e) and inter._energy > 10*GeV:
							high_e = True
							geo.scatter3d(self._tree.voxx, self._tree.voxy, self._tree.voxz, self._tree.voxe, "event_" + str(k), inter.particle_string() + " KE=" + str(round(inter._energy/GeV) ), False)
						if low_e and high_e:
							break

	def stack(self):
		flavor_stack = []
		current_stack = []
		class_stack = []

		ENERGY_NBINS = 500
		ENERGY_MIN = 0.00 
		ENERGY_MAX = 20.0

		total_energy = root.TH1D("sum", "sum", ENERGY_NBINS, ENERGY_MIN, ENERGY_MAX)
		total_energy.SetLineColor(1)
		total_energy.Sumw2(True)

		for k, flavor in enumerate(self.flavors):
			flavor_stack.append(root.TH1D(flavor, flavor, ENERGY_NBINS, ENERGY_MIN, ENERGY_MAX))
		for k, current in enumerate(self.currents):
			current_stack.append(root.TH1D(current, current, ENERGY_NBINS, ENERGY_MIN, ENERGY_MAX))
		for k, _class in enumerate(self.classes):
			class_stack.append(root.TH1D(_class, _class, ENERGY_NBINS, ENERGY_MIN, ENERGY_MAX))

		for k in range(self._tree.GetEntries()):
			self._tree.GetEntry(k)
			inter = physics.interaction(self._tree.fsPDG, self._tree.fsKE, self._tree.voxe)

			for k, flavor in enumerate(self.flavors):
				if inter._flavor == flavor:
					flavor_stack[k].Fill(inter._energy/GeV)

			for k, current in enumerate(self.currents):
				if inter._current == current:
					current_stack[k].Fill(inter._energy/GeV)

			for k, _class in enumerate(self.classes):
				if inter._class == _class:
					class_stack[k].Fill(inter._energy/GeV)

			total_energy.Fill(inter._total_dep/GeV)

		
		root.gStyle.SetOptStat(1)

		fstack = root.THStack("flavor", "Energy by Neutrino Flavor")
		cstack = root.THStack("current", "Energy by W/Z Current")
		tstack = root.THStack("class", "Energy Deposited by Interaction Type")

		flegend = root.TLegend(0.60, 0.65, 0.90, 0.90)
		clegend = root.TLegend(0.60, 0.65, 0.90, 0.90)
		tlegend = root.TLegend(0.60, 0.65, 0.90, 0.90)

		integral = 0.0
		for k, plot in enumerate(flavor_stack):
			integral = integral + plot.Integral()
		for k, plot in enumerate(flavor_stack):
			plot.SetFillColor(5*k+52)
			plot.SetLineColor(5*k+52)
			fstack.Add(plot)
			flegend.AddEntry(plot, self.flavors[k] + " " + str(round(100.0*plot.Integral()/integral, 2)) + "%")

		for k, plot in enumerate(current_stack):
			plot.SetFillColor(5*k+52)
			plot.SetLineColor(5*k+52)
			clegend.AddEntry(plot, self.currents[k] + " " + str(round(100.0*plot.Integral()/integral, 2)) + "%")
			cstack.Add(plot)
		for k, plot in enumerate(class_stack):
			plot.SetFillColor(3*k+51)
			plot.SetLineColor(3*k+51)
			tlegend.AddEntry(plot, self.classes[k] + " " + str(round(100.0*plot.Integral()/integral, 2)) + "%")
			tstack.Add(plot)


		root.gStyle.SetOptStat(2)

		c1 = root.TCanvas("c1", "c1")
		total_energy.Draw()
		c1.Print("total.pdf", ".pdf")

		fcanv = root.TCanvas("cf", "cf")
		fcanv.cd()
		fcanv.SetLogy()
		fstack.Draw()
		fstack.GetXaxis().SetTitle("Total Final State Energy [GeV]")
		fstack.GetXaxis().CenterTitle()
		flegend.Draw("SAME")
		fcanv.Print("flavor_stack.pdf", ".pdf")

		ccanv = root.TCanvas("cc", "cc")
		ccanv.cd()
		ccanv.SetLogy()
		cstack.Draw()
		cstack.GetXaxis().SetTitle("Total Final State Energy [GeV]")
		cstack.GetXaxis().CenterTitle()
		clegend.Draw("SAME")
		ccanv.Print("current_stack.pdf", ".pdf")

		tcanv = root.TCanvas("ct", "ct")
		tcanv.cd()
		tcanv.SetLogx()
		tcanv.SetLogy()


		tstack.Draw()
		tstack.GetXaxis().SetTitle("Total Energy Deposit [GeV]")
		tstack.GetXaxis().CenterTitle()
		tlegend.Draw("SAME")
		tcanv.Print("type_stack.pdf", ".pdf")

	











			



	


