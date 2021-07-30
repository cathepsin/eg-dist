import Sequence
import Chain
import ResidueCentroid
import HeptadRepeat
import FPTools as fpt
import tkinter as tk
import sys
from tkinter.filedialog import askopenfilename

root = tk.Tk()
root.withdraw()
f_pdb = askopenfilename(title="Select a .pdb file", filetypes=[('Protein Database File', '*.pdb')
    , ('Protein Database File (old)', '*.ent')])

try:
    f_p = open(f_pdb)
except FileNotFoundError:
    sys.exit("Must select a file")

protSeq = Sequence.ProteinSequence()
protSeq.parsePDB(f_p)
seqChains = Chain.Chain(protSeq)
centroid = ResidueCentroid.CentroidFinder()
rCentroids = []
for aa in protSeq.sequence:
    rCentroids.append(centroid.GetCentroid(aa))

pmt = "Select associated .socket file for " + fpt.CutPath(f_pdb)
f_sock = askopenfilename(title=pmt, filetypes=[('Socket File', '*.socket')])
try:
    f_s = open(f_sock)
except FileNotFoundError:
    sys.exit("Must select a file")

heptad = HeptadRepeat.Heptad()
heptad.ParseSocket(f_s)

# print(protSeq.GetAtomDist(protSeq.sequence[14].atoms[1], protSeq.sequence[128].atoms[1]))
print("Done")
