from prody import *
from pylab import *

import sys

ion()

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

filepath = input("Input file path to .pdb file ")
try:
    f = open(filepath)
except (FileNotFoundError, OSError):
    filepath = fStobS(filepath)
    try:
        f = open(filepath)
    except (FileNotFoundError, OSError):
        sys.stderr.write("File not found. Make sure your file is at the correct location and is a .pdb\nAborted")
        sys.exit(1) #Abort



protein = prody.parsePDBStream(f)
prody.showProtein(protein)
print("wow")

