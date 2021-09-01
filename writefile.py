import ntpath
import os

import EGDist
import HeptadRepeat


class FileWriter:
    def __init__(self, *argv):
        self.arg = argv
        self.summary = ""
        print("Putting master file in: ", os.getcwd())
        self.outfile = open("Master_Summary.csv", "w")

    def SetWriteFile(self, file):
        self.writefile = file
        self.writefile.write("Protein from " + ntpath.basename(file.name) + ",\n")
        self.writefile.write("Pairs:,CA,CB,Centroid,\n")

    def SetArgs(self, *argv):
        self.arg = argv

    #Write all data to a file from given argv. Fill summary
    def writeFile(self):
        egdist = EGDist.EGDist()
        heptad = HeptadRepeat.Heptad()
        errorMessage = ""
        for clss in self.arg:
            errorMessage += clss.message
            if isinstance(clss, type(egdist)):
                self.writefile.write("," + str(clss.distances["CA"]) + "," + str(clss.distances["CB"]) + "," + str(clss.distances["Centroid"]) + ",")
                self.summary += self.writefile.name + "," + str(clss.distances["CA"]) + "," + str(clss.distances["CB"]) + "," + str(clss.distances["Centroid"]) + ","
            if isinstance(clss, type(heptad)):
                oligomeric_state = clss.oligomeric_state
        self.writefile.write(oligomeric_state + ",\n")
        self.summary += oligomeric_state + ",\n"
        self.writefile.write("\nWarnings:\n" + errorMessage)

    #Write to master csv file
    def WriteMaster(self, num):
        self.outfile.write(f"Master Summary for {num} proteins\n")
        self.outfile.write("File name,CA e-g distance,CB e-g distance,Centroid e-g distance,Oligomeric State\n")
        self.outfile.write(self.summary)
