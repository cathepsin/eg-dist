# Converts backslashes to forward slashes. Allows filepath lookup in linux and windows as needed
def bStofS(fp):
    retFP = ""
    for i in range(len(fp)):
        if fp[i] == "\\":
            retFP += "/"
        else:
            if (i == 0 or i == len(fp) - 1) and fp[i] == "\"":
                continue
            retFP += fp[i]

    return retFP


# Cuts the file path from a file name
def CutPath(str):
    retStr = str
    while retStr.find('/') != -1:
        retStr = retStr[retStr.find('/') + 1:]
    return retStr