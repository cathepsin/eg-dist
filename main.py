import Sequence
import Chain
import ResidueCentroid
import HeptadRepeat
import EGDist
import FPTools as fpt
import tkinter as tk
import sys
from tkinter.filedialog import askopenfilename

root = tk.Tk()
root.withdraw()
f_pdb = askopenfilename(title="Select a .pdb file", filetypes=[
    ('Molecular Topology', '*.mmol')
    , ('Protein Database File (old) (Not yet supported)', '*.ent')
    , ( 'Protein Database File (Not yet supported)', '*.pdb')])

try:
    f_p = open(f_pdb)
except FileNotFoundError:
    sys.exit("Must select a file")

protSeq = Sequence.ProteinSequence()
protSeq.parsePDB(f_p)
seqChains = Chain.Chain(protSeq)
f_p.seek(0)
seqChains.Symmetry(f_p, protSeq)
for key in seqChains.matrix:
    protSeq.GeneratePair(seqChains.copyChains, seqChains.matrix[key])
centroid = ResidueCentroid.CentroidFinder()
rCentroids = []
for aa in protSeq.sequence:
    rCentroids.append(centroid.GetCentroid(aa))

pmt = "Select associated .socket file for " + fpt.CutPath(f_pdb)
f_sock = askopenfilename(title=pmt, filetypes=[('Socket File', '*.short.socket')])
try:
    f_s = open(f_sock)
except FileNotFoundError:
    sys.exit("Must select a file")
heptad = HeptadRepeat.Heptad()
heptad.ParseSocket(f_s)
#TODO Take centroid of only the helices. No need to calculate centroid of other residues.
distances = EGDist.EGDist()
distances.GetDistances(protSeq.sequence, seqChains, rCentroids, heptad)
# print(protSeq.GetAtomDist(protSeq.sequence[14].atoms[1], protSeq.sequence[128].atoms[1]))
print("Done")
