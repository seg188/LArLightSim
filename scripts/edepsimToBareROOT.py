#!/usr/bin/env python

import os
import ROOT
from optparse import OptionParser
import subprocess
from array import array

MAXHITS = 1000000
MAXPARTICLES = 10000
MAXINT = 100
t_ev = array( 'i', [0] )
t_vtxx = array( 'd', [0.] )
t_vtxy = array( 'd', [0.] )
t_vtxz = array( 'd', [0.] )
t_nvox = array( 'i', [0] )
t_voxx = array( 'd', [0.]*MAXHITS )
t_voxy = array( 'd', [0.]*MAXHITS )
t_voxz = array( 'd', [0.]*MAXHITS )
t_voxe = array( 'd', [0.]*MAXHITS )
t_vox_dedx = array( 'd', [0.]*MAXHITS )
t_voxTID = array( 'd', [0]*MAXHITS )
t_intID = array( 'd', [0]*MAXHITS )
t_nFSP = array( 'i', [0] )
t_fsPDG = array( 'i', [0]*MAXPARTICLES )
t_fsKE = array( 'd', [0.]*MAXPARTICLES )

offset = ROOT.TVector3( 0., 73.94, 479.25 )
dimensions = ROOT.TVector3( 713.1, 300., 507.3 ) # cm
padpitch = 0.4 # cm
timeres = 0.4 # cm, resolution is better in x dimension because timing

def getVoxyl( pos ):
    x = int( (pos.x() + dimensions.x()/2.) / timeres )
    y = int( (pos.y() + dimensions.y()/2.) / padpitch )
    z = int( pos.z() / padpitch )
    return int(x), int(y), int(z)

def loop( events, tgeo, hdf5file ):

    event = ROOT.TG4Event()
    events.SetBranchAddress("Event",ROOT.AddressOf(event))

    N = events.GetEntries()
    
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
            vtx = 0.1*(vertex.Position.Vect())
            vtx += offset
            
            t_vtxx[0] = vtx.x()
            t_vtxy[0] = vtx.y()
            t_vtxz[0] = vtx.z()
            # We probably want to save information about the true neutrino interaction, so that we can 
            # see how the reco does for different types of interactions
            t_nFSP[0] = 0
            for ipart,particle in enumerate(vertex.Particles):
                mom = particle.Momentum
                pdg = particle.PDGCode
                ke = mom.E() - mom.M()
                if pdg >= 9999:
                    # don't count nuclear fragments, bindinos, etc.
                    continue
                t_fsPDG[t_nFSP[0]] = pdg
                t_fsKE[t_nFSP[0]] = ke
                t_nFSP[0] += 1

            ArCube_hits = []
            for det in event.SegmentDetectors:
                if det.first == "ArgonCube":
                    ArCube_hits = det.second

            voxyl_energy = {}
            for k, edep in enumerate(ArCube_hits):
                if edep.EnergyDeposit < 0.01: # Some energy threshold, this is 10 keV
                    continue

                node = tgeo.FindNode( edep.Start.X(), edep.Start.Y(), edep.Start.Z())
                if "volLArActive" not in node.GetName():
                    continue

                hStart = ROOT.TVector3( edep.Start.X()/10., edep.Start.Y()/10., edep.Start.Z()/10. )
                hStop = ROOT.TVector3( edep.Stop.X()/10., edep.Stop.Y()/10., edep.Stop.Z()/10. )
                hStart += offset # (0,0,0) is in the middle of the upstream face of the active volume
                hStop += offset

                hDir = (hStop-hStart).Unit()
                mag = (hStop-hStart).Mag()
                # Chop the hit up
                n_chop = int(mag/0.01) + 1 # round up
                step_length = mag/n_chop
                dedx = edep.EnergyDeposit/step_length

                for i in range(n_chop):
                    hPos = hStart + ((i*step_length)*hDir)
                    x,y,z = getVoxyl( hPos ) # x,y,z are integer coordinates of voxyls
                    if (x,y,z) not in voxyl_energy:
                        voxyl_energy[(x,y,z)] = edep.EnergyDeposit/n_chop
                    else:
                        voxyl_energy[(x,y,z)] += edep.EnergyDeposit/n_chop


            t_nvox[0] = 0
            for voxyl in voxyl_energy:
                t_voxx[t_nvox[0]] = -dimensions.x()/2. + voxyl[0]*timeres
                t_voxy[t_nvox[0]] = -dimensions.y()/2. + voxyl[1]*padpitch
                t_voxz[t_nvox[0]] = voxyl[2]*padpitch
                t_voxe[t_nvox[0]] = voxyl_energy[voxyl]
                t_vox_dedx[t_nvox[0]] = dedx
                t_nvox[0] += 1

            tree.Fill()

if __name__ == "__main__":

    ROOT.gROOT.SetBatch(1)

    parser = OptionParser()
    parser.add_option('--outfile', help='Output file name', default="out.root")
    parser.add_option('--topdir', help='Input file top directory', default="~/hex/berk/edep/neutrino/DetEnclosure")
    parser.add_option('--first_run', type=int, help='First run number', default=1001)
    parser.add_option('--last_run', type=int, help='Last run number', default=1001)
    parser.add_option('--rhc', action='store_true', help='Reverse horn current', default=False)
    parser.add_option('--geom',help='top volume of interactions', default="DetEnclosure")
    parser.add_option('--grid',action='store_true', help='Grid mode')

    (args, dummy) = parser.parse_args()

    tgeo = None

    neutrino = "neutrino" if not args.rhc else "antineutrino"
    horn = "FHC" if not args.rhc else "RHC"

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
    tree.Branch( "vox_dedx", t_vox_dedx, "voxe[nvox]/D" )
    tree.Branch( "nFSP", t_nFSP, "nFSP/I" )
    tree.Branch( "fsPDG", t_fsPDG, "fsPDG[nFSP]/I" )
    tree.Branch( "fsKE", t_fsKE, "fsKE[nFSP]/D" )

    for run in range( args.first_run, args.last_run+1 ):
        fname = None
        if args.grid:
            fname = "%s.%d.edepsim.root" % (neutrino, run)
        else:
            fname = "%s/EDep/%s/%s/%s.%d.edepsim.root" % (args.topdir, horn, args.geom, neutrino, run)
        if not os.access( fname, os.R_OK ):
            print ("Can't access file: %s" % fname)
            #continue
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


