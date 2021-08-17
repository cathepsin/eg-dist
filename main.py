import Sequence
import Chain
import ResidueCentroid
import HeptadRepeat
import EGDist
import filemanager
import CustomExceptions
import FPTools as fpt
import tkinter as tk
import sys
import os
from tkinter.filedialog import askopenfilename
import ntpath
import signal

ACCEPTED_EXTENSIONS = [".mmol", ".ent", "pdb", ".short.socket"]
def CheckExtensions(val):
    for ext in ACCEPTED_EXTENSIONS:
        loc = val.find(ext)
        if val[loc:] == ext:
            return True
    return False





currdir = os.getcwd()
def signal_handler(signal, frame):
    os.chdir(currdir)
    print(os.getcwd())
    print("here")
    sys.exit(-1)
signal.signal(signal.SIGINT, signal_handler)


ACCEPTED_FILES = [
    ('Molecular Topology', '*.mmol')
    , ('Protein Database File (old) (Not yet supported)', '*.ent')
    , ( 'Protein Database File (Not yet supported)', '*.pdb')
]
ACCEPTED_SOCKET = [('Socket File', '*.short.socket')]

ERROR_MESSAGE = "Must select a file"
PDB_PROMPT = "Select a .pdb file"

try:
    if len(sys.argv) == 1:
        print("tkinter")
        root = tk.Tk()
        root.withdraw()
        f_pdb = askopenfilename(title=PDB_PROMPT, filetypes=ACCEPTED_FILES)
        pmt = "Select associated .socket file for " + fpt.CutPath(f_pdb)
        f_sock = askopenfilename(title=pmt, filetypes=ACCEPTED_SOCKET)
        filelist = [filemanager.TwoFiles(f_pdb, f_sock)]
    elif sys.argv[1] == "-file":
        print("give me two files!")
        if CheckExtensions(sys.argv[2]) and CheckExtensions(sys.argv[3]):
            if sys.argv[2].endswith(".short.socket") or sys.argv[3].endswith(".short.socket") and not (sys.argv[2].endswith(".short.socket") and sys.argv[3].endswith(".short.socket")):
                print(sys.argv[2], sys.argv[3])
                filelist = [filemanager.TwoFiles(sys.argv[2], sys.argv[3])]
    elif sys.argv[1] == "-dir":
        if not os.path.isdir(sys.argv[2]):
            raise CustomExceptions.NotDirectory
        print("give me a directory")
        filemanager.GetFileList(sys.argv[2])
        filemanager.Organize()
        filelist = [file for file in os.listdir() if os.path.isdir(file)]

    for file in filelist:
        print(os.path.join(os.getcwd(), file))
        os.chdir(os.path.join(os.getcwd(), file))
        for ptr in os.listdir():
            if ptr.find(".socket") == -1:
                f_pdb = ptr
            else:
                f_sock = ptr

    #f_pdb = askopenfilename(title=PDB_PROMPT, filetypes=ACCEPTED_FILES)

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
        #f_sock = askopenfilename(title=pmt, filetypes=ACCEPTED_SOCKET)
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

        #TODO get oligomeric state for csv
        print(ntpath.basename(f_p.name))
        outfile = open(ntpath.basename(f_p.name) + ".csv", "w")
        outfile.write("Protein from " + ntpath.basename(f_p.name) + ",\n")
        outfile.write("Pairs:,CA,CB,Centroid,\n")
        print("Done")
        outfile.close()
        os.chdir("..")
finally:
    os.chdir(currdir)