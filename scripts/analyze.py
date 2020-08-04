import physics
import read_data as rd 
import ROOT as root

data_file = root.TFile.Open("~/hex/berk/nersc/data/stephenBig.root")
tree = data_file.Get("tree")

tree.GetEntry(398)
inter = physics.interaction(tree.fsPDG, tree.fsKE)
sup = rd.supervisor()
d = rd.data(tree, sup)
d.draw_characteristic_events()

print(inter._particles)