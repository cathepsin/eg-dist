import tkinter as tk
import Sequence
import Chain
from tkinter.filedialog import askopenfilename


def fStobS(fp):
    retFP = ""
    for i in range(len(fp)):
        if fp[i] == "\\":
            retFP += "/"
        else:
            if (i == 0 or i == len(fp) - 1) and fp[i] == "\"":
                continue
            retFP += fp[i]

    return retFP


def CutPath(str):
    # Cuts the file path from a file name
    retStr = str
    while retStr.find('/') != -1:
        retStr = retStr[retStr.find('/') + 1:]
    return retStr


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


pmt = "Select associated .socket file for " + CutPath(f_pdb)
f_sock = askopenfilename(title=pmt, filetypes=[('Socket File', '*.socket')])
try:
    f_s = open(f_sock)
except FileNotFoundError:
    sys.exit("Must select a file")

# print(protSeq.GetAtomDist(protSeq.sequence[14].atoms[1], protSeq.sequence[128].atoms[1]))
print("Done")
