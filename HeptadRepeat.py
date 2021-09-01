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
        self.message = ""
        self.oligomeric_state = "Unknown"
        self.heptadInfo = []

    #Parse a .short.socket file to extract information about a heptad repeat
    def ParseSocket(self, file):
        rFile = file.readlines()
        i = 0
        while i < len(rFile):
            #Get oligomeric state from .short.socket file. Gets only first instance of described oligomeric state
            if rFile[i].find("stranded") != -1 and self.oligomeric_state == "Unknown":
                for word in rFile[i].split():
                    if word.find("stranded") != -1:
                        self.oligomeric_state = word
                        for i in range(len(self.oligomeric_state)):
                            if not self.oligomeric_state[i].isalnum():
                                self.oligomeric_state = self.oligomeric_state[:i] + " " + self.oligomeric_state[i + 1:]
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

    #Get the entire paragraph containing the heptad repeat information for easier handling
    def ParseParagraph(self, para):
        lines = para.splitlines()
        lines[0] = lines[0].replace('-',' ')
        lines[0] = lines[0].replace(':',' ')
        firstline = lines[0].split()
        chain = firstline[len(firstline) - 1]
        try:
            range = [int(firstline[len(firstline) - 3]), int(firstline[len(firstline) - 2])]
        except ValueError:
            range = [int(firstline[len(firstline) - 2]), int(firstline[len(firstline) - 1])]
        sequence = lines[2].split()[1]
        lines[3] = lines[3].replace(' ', '-').replace('-',' ',1)
        register = lines[3].split()[1]
        if chain == str(range[1]):#If socket file was produced from pdf file with missing chain labels
            chain = "$"
        return chain, range, sequence, register

    #AminoAcids getter
    def GetList(self):
        return self.AAs

#Convert one-letter code to three-letter code. Note this is defined outside of the class
def OneToThree(AA):
    object = Heptad()
    return object.GetList()[AA]