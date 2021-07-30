class Heptad:
    def __init__(self):
        self.AAs = {
            "R": 'ARG',
            "H": 'HIS',
            "K": 'LYS',
            "D": 'ASP',
            "E": 'GLU',
            "S": 'SER',
            "T": 'THR',
            "N": 'ASN',
            "Q": 'GLN',
            "C": 'CYS',
            "U": 'SEC',
            "G": 'GLY',
            "P": 'PRO',
            "A": 'ALA',
            "V": 'VAL',
            "I": 'ILE',
            "L": 'LEU',
            "M": 'MET',
            "F": 'PHE',
            "Y": 'TYR',
            "W": 'TRP',
            "O": 'PYL'
        }

    def ParseSocket(self, file):
        rFile = file.readlines()
        i = 0
        while i < len(rFile):
            if rFile[i].find("assigning heptad to helix") == 0:
                paragraph = ""
                while rFile[i] != "\n":
                    paragraph += rFile[i]
                    if rFile[i].strip() == "Finished":
                        break
                    i += 1

                print(paragraph)

            if rFile[i].strip() == "Finished":
                break
            i += 1

