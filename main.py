import Sequence
import Chain
import ResidueCentroid
import HeptadRepeat
import EGDist
import FPTools as fpt
import tkinter as tk
import sys
from tkinter.filedialog import askopenfilename

ACCEPTED_FILES = [
    ('Molecular Topology', '*.mmol')
    , ('Protein Database File (old) (Not yet supported)', '*.ent')
    , ( 'Protein Database File (Not yet supported)', '*.pdb')
]
ACCEPTED_SOCKET = [('Socket File', '*.short.socket')]
ERROR_MESSAGE = "Must select a file"
PDB_PROMPT = "Select a .pdb file"

root = tk.Tk()
root.withdraw()
f_pdb = askopenfilename(title=PDB_PROMPT, filetypes=ACCEPTED_FILES)

try:
    f_p = open(f_pdb)
except FileNotFoundError:
    sys.exit(ERROR_MESSAGE)

protSeq = Sequence.ProteinSequence()
protSeq.parsePDB(f_p)
seqChains = Chain.Chain(protSeq)
f_p.seek(0)
seqChains.Symmetry(f_p, protSeq)
for key in seqChains.matrix:
    protSeq.GeneratePair(seqChains.copyChains, seqChains.matrix[key])
pmt = "Select associated .socket file for " + fpt.CutPath(f_pdb)
f_sock = askopenfilename(title=pmt, filetypes=ACCEPTED_SOCKET)
try:
    f_s = open(f_sock)
except FileNotFoundError:
    sys.exit(ERROR_MESSAGE)
heptad = HeptadRepeat.Heptad()
heptad.ParseSocket(f_s)
distances = EGDist.EGDist()
distances.GetRangeSequence(heptad.heptadInfo, seqChains)

centroid = ResidueCentroid.CentroidFinder()
rCentroids = []
for helix in distances.helices:
    for aa in helix:
        aa.SetCentroid(centroid.GetCentroid(aa))

distances.GetDistances()
print("Done")
