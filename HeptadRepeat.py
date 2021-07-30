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

        self.heptadInfo = []

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

                self.heptadInfo.append(self.ParseParagraph(paragraph))

            if rFile[i].strip() == "Finished":
                break
            i += 1


    def ParseParagraph(self, para):
        lines = para.splitlines()
        lines[0] = lines[0].replace('-',' ')
        lines[0] = lines[0].replace(':',' ')
        firstline = lines[0].split()
        chain = firstline[len(firstline) - 1]
        range = [firstline[len(firstline) - 3], firstline[len(firstline) - 2]]
        sequence = lines[2].split()[1]
        lines[3] = lines[3].replace(' ', '-').replace('-',' ',1)
        register = lines[3].split()[1]
        return chain, range, sequence, register

    def GetList(self):
        return self.AAs

def OneToThree(AA):
    object = Heptad()
    return object.GetList()[AA]