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