#!/usr/bin/env python

import os
import ROOT
from optparse import OptionParser
import subprocess
from array import array
import numpy as np
E_NBINS = 1000
E_MIN = 0.001
E_MAX = 100.0
P_NBINS = 1000
P_MIN = 0.10
P_MAX = 5000
R_NBINS = 1000
R_MIN = 0.99
R_MAX = 1.00
mu_mass = 104.0
calc_eloss_1d =  ROOT.TH1D("eloss_1d", "dedx", E_NBINS, E_MIN, E_MAX)
plr = ROOT.TH1D("secondary_ratio", "secondary deposit ratio",R_NBINS, R_MIN, R_MAX )
gr_steps = ROOT.TH1D("steps2", "step length process 2", E_NBINS, E_MIN, E_MAX)
gr_steps4 = ROOT.TH1D("steps4", "step length process 2", E_NBINS, E_MIN, E_MAX)
gr_dedx = ROOT.TH1D("dedx", "dedx", E_NBINS, E_MIN, E_MAX)
gr_dedx_primary = ROOT.TH1D("dedx_prim", "dedx", E_NBINS, E_MIN, E_MAX)
gr_dedx_secondary = ROOT.TH1D("dedx_sec", "dedx", E_NBINS, E_MIN, E_MAX)
gr_eloss = ROOT.TH2D("eloss", "de/dx vs. p for incident muon", P_NBINS, P_MIN, P_MAX, E_NBINS, E_MIN, E_MAX)
gr_eloss_dep = ROOT.TH2D("eloss_dep", "de/dx vs. p for incident muon", P_NBINS, P_MIN, P_MAX, E_NBINS, E_MIN, E_MAX)

MAXHITS = 1000000
MAXPARTICLES = 10000
MAXINT = 100
output_es = [1001]#[100, 300, 500, 750, 1000, 1500, 2000, 1250, 1750, 200, 50, 1100, 20, 900, 1400, 1600, 1850, 125]
gr_z = ROOT.TGraph(MAXHITS)
t_ev = array( 'i', [0] )
t_vtxx = array( 'd', [0.] )
t_vtxy = array( 'd', [0.] )
t_vtxz = array( 'd', [0.] )
t_nvox = array( 'i', [0] )
t_voxx = array( 'd', [0.]*MAXHITS )
t_voxy = array( 'd', [0.]*MAXHITS )
t_voxz = array( 'd', [0.]*MAXHITS )
t_voxe = array( 'd', [0.]*MAXHITS )
t_dedx = array( 'd', [0.]*MAXHITS )
t_voxTID = array( 'd', [0]*MAXHITS )
t_intID = array( 'd', [0]*MAXHITS )
t_nFSP = array( 'i', [0] )
t_fsPDG = array( 'i', [0]*MAXPARTICLES )
t_fsKE = array( 'd', [0.]*MAXPARTICLES )

proccesses = [91, 2, 401]
proc_plots = []
for k in range(len(proccesses)):
    proc_plots.append(ROOT.TH1D(str(proccesses[k]), str(proccesses[k]) , E_NBINS, E_MIN, E_MAX ))


offset = ROOT.TVector3( 0., 73.94, 479.25 )
dimensions = ROOT.TVector3( 713.1, 300., 507.3 ) # cm
padpitch = 0.4 # cm
timeres = 0.4 # cm, resolution is better in x dimension because timing

def getVoxyl( pos ):
    x = int( (pos.x() + dimensions.x()/2.) / timeres )
    y = int( (pos.y() + dimensions.y()/2.) / padpitch )
    z = int( pos.z() / padpitch )
    return int(x), int(y), int(z)
def dist(tlorentz):
    x, y, z = tlorentz.X(), tlorentz.Y(), tlorentz.Z()
    return np.sqrt(x*x + y*y + z*z)

def loop( events, tgeo, hdf5file ):

    event = ROOT.TG4Event()
    events.SetBranchAddress("Event",ROOT.AddressOf(event))

    N = events.GetEntries()
    print(N)
    event.Print()
    for ient in range(N):
        events.GetEntry(ient) # Load edep-sim event

        info = ROOT.ProcInfo_t()
        ROOT.gSystem.GetProcInfo(info)
        vm = info.fMemVirtual
        rm = info.fMemResident

       # if ient % 100 == 0:
     #       print "Event %d of %d, VM = %1.2f MB, RM = %1.2f MB" % (ient,N,vm/1000.,rm/1000.)

        t_ev[0] = ient
        # There is only ever one primary vertex per event, but in principle there could be more, so we have to do this
        for ivtx,vertex in enumerate(event.Primaries):

            # Save the vertex position so we can determine if it's in the fiducial volume or not
            vtx = 0.1*(vertex.GetPosition().Vect())
            vtx += offset
            
            t_vtxx[0] = vtx.x()
            t_vtxy[0] = vtx.y()
            t_vtxz[0] = vtx.z()

            # We probably want to save information about the true neutrino interaction, so that we can 
            # see how the reco does for different types of interactions
            t_nFSP[0] = 0
            for ipart,particle in enumerate(vertex.Particles):
                mom = particle.GetMomentum()
                pdg = particle.GetPDGCode()
                ke = mom.E() - mom.M()
                if pdg >= 9999:
                    # don't count nuclear fragments, bindinos, etc.
                    continue
                t_fsPDG[t_nFSP[0]] = pdg
                t_fsKE[t_nFSP[0]] = ke
                t_nFSP[0] += 1
        graphed_z = False

        for itr, traj in enumerate(event.Trajectories):
            pdg = traj.GetPDGCode()
            if pdg == 13 or pdg == -13:
                L = len(traj.Points)

                for ipt, pt in enumerate(traj.Points):
                    #LOOK THROUGH TRAJECTORY POINTS HERE
                    if ipt == 0:
                        continue


                    pt_fwd = pt.GetPosition()
                    pt_back = traj.Points[ipt-1].GetPosition()

                    x = pt_fwd.X() - pt_back.X()
                    y = pt_fwd.Y() - pt_back.Y()
                    z = pt_fwd.Z() - pt_back.Z()


                    step_length = np.sqrt(x*x + y*y + z*z)*0.10
                    E_i = 0.
                    E_f = 0.
                    subpc = (pt.GetProcess())
                    if step_length > 0:

                        mom = pt.GetMomentum().Mag()
                        mom_i = traj.Points[ipt-1].GetMomentum().Mag()

                        E_i = np.sqrt(mom_i*mom_i + mu_mass*mu_mass)
                        E_f = np.sqrt(mom*mom + mu_mass*mu_mass)

                        if subpc == 401 or step_length < 10:
                            gr_eloss.Fill(mom, (E_i - E_f)/step_length)
                            gr_steps4.Fill(step_length)

                        if subpc == 2:
                            gr_steps.Fill(step_length)

                        prc = (pt.GetProcess())
                        subpc = (pt.GetSubprocess())
                     
                        for ipr, proc in enumerate(proccesses):
                            if subpc == proc:

                                proc_plots[ipr].Fill((E_i - E_f)/step_length)
                        if prc == 2 or prc == 7:
                     
                            calc_eloss_1d.Fill( 1000.0*(E_i - E_f)/step_length )

                if L > 100:
                    graphe_z =True
                            



        ArCube_hits = []
        for n, det in enumerate(event.SegmentDetectors):
            #if det.first == "ArgonCube":
        


            if n == 0:
                ArCube_hits = det.second
                #print(det.second)

        for k, edep in enumerate(ArCube_hits):
            _id = (edep.GetPrimaryId())
            if _id == 0:
                d = dist(edep.GetStart())

                
                gr_dedx.Fill(10.0 * edep.GetEnergyDeposit() / edep.GetTrackLength())
                
                gr_dedx_primary.Fill(10.0 * edep.GetEnergyDeposit() / edep.GetTrackLength())
                gr_dedx_secondary.Fill(10.0 * edep.GetSecondaryDeposit() / edep.GetTrackLength())

            if edep.GetEnergyDeposit() > 0:
                plr.Fill(edep.GetSecondaryDeposit()/edep.GetEnergyDeposit())
             
              
        tree.Fill()

if __name__ == "__main__":


    ROOT.gROOT.SetBatch(1)

    parser = OptionParser()
    parser.add_option('--outfile', help='Output file name', default="out.root")
    parser.add_option('--file', help='edep-sim Data File', default="/home/stephen/hex/berk/edep/new/edep-sim/output/")

    (args, dummy) = parser.parse_args()

    tgeo = None

    #SETTING OUTPUT FILE AND TREE BRANCHES

    outfile = ROOT.TFile( args.outfile, "RECREATE" )
    tree = ROOT.TTree( "tree", "tree" )
    tree.Branch( "ev", t_ev, "ev/I" )
    tree.Branch( "vtxx", t_vtxx, "vtxx/D" )
    tree.Branch( "vtxy", t_vtxy, "vtxy/D" )
    tree.Branch( "vtxz", t_vtxz, "vtxz/D" )
    tree.Branch( "nvox", t_nvox, "nvox/I" )
    tree.Branch( "voxx", t_voxx, "voxx[nvox]/D" )
    tree.Branch( "voxy", t_voxy, "voxy[nvox]/D" )
    tree.Branch( "voxz", t_voxz, "voxz[nvox]/D" )
    tree.Branch( "voxe", t_voxe, "voxe[nvox]/D" )
    tree.Branch( "dedx", t_dedx, "dedx[nvox]/D" )
    tree.Branch( "nFSP", t_nFSP, "nFSP/I" )
    tree.Branch( "fsPDG", t_fsPDG, "fsPDG[nFSP]/I" )
    tree.Branch( "fsKE", t_fsKE, "fsKE[nFSP]/D" )

 
    folder_name = args.file

    for e in output_es:
        fname = folder_name + "output" + str(int(e)) + ".root"

        tf = ROOT.TFile( fname )
        events = tf.Get( "EDepSimEvents" )

        if tgeo is None:
            tf.MakeProject("EDepSimEvents","*","RECREATE++")
            tgeo = tf.Get( "EDepSimGeometry" )

        print ("Looping over: %s" % fname)
        # Loop over one edep-sim input file
        outfile.cd()
        loop( events, tgeo, tree )
        tf.Close()

    outfile.cd()
    tree.Write()
    gr_dedx.Write() 
    gr_eloss.Write()
    gr_dedx_secondary.Write()
    gr_dedx_primary.Write()
    gr_steps.Write()
    gr_steps4.Write()
    gr_z.Write()
    plr.Write()
    calc_eloss_1d.Write()

    for k, plot in enumerate(proc_plots):
        plot.Write()


