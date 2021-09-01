import math
import numpy as np

class EGDist:
    def __init__(self):
        self.message = ""

    class AAPair: #Class to store data about an e-g pair
        def __init__(self, e, g):
            self.e = e
            self.g = g
        def setVD(self, val):
            self.zdist = val

        def setEPoint(self, pt):
            self.EPoint = pt

        def setGPoint(self, pt):
            self.GPoint = pt

        #Get distance from two points
        def GetDist(self, a1, a2):
            return ((a1[0] - a2[0]) ** 2 + (a1[1] - a2[1]) ** 2 + (a1[2] - a2[2]) ** 2) ** (1 / 2)

        #Get distance between alpha carbons
        def CADist(self):
            # Get CAs
            ca1 = None
            ca2 = None
            for atom in self.e.atoms:
                if atom.id == 'CA':
                    ca1 = atom
                    break
            for atom in self.g.atoms:
                if atom.id == 'CA':
                    ca2 = atom
                    break
            return self.GetDist(ca1.location, ca2.location)

        #Get distance between beta carbons
        def CBDist(self):
            # Get CBs
            cb1 = None
            cb2 = None
            for atom in self.e.atoms:
                if atom.id == 'CB':
                    cb1 = atom
                    break
            for atom in self.g.atoms:
                if atom.id == 'CB':
                    cb2 = atom
                    break
            return self.GetDist(cb1.location, cb2.location)

        #Get centroid distance
        def CentDist(self):
            return self.GetDist(self.e.centroid, self.g.centroid)

        # Class dunder methods
        def __eq__(self, other):
            if isinstance(self, type(other)):
                return self.e == other.e and self.g == other.g
            else:
                return False

        def __ne__(self, other):
            if not isinstance(self, type(other)):
                return True
            elif self.e != other.e or self.g != other.g:
                return True
            return False

        def __repr__(self):
            return f"{self.e} || {self.g}"

        def __str__(self):
            return f"{self.e} || {self.g}"

    #Get distances for all potential pairings. Keep only likely pairings
    def GetDistances(self):
        self.pairs = list()
        for helix in self.helices:
            #Select a helix. Compare it only to other helices
            exclusionList = [hel for hel in self.helices if hel != helix]
            for residue in helix:
                pair = self.GetPair(residue, exclusionList)
                #If a pair is incorrectly created with a missing residue, continue to the next possible pairing
                if len(pair) == 1:
                    continue
                if pair[0] == None:
                    continue
                if pair[1] == None:
                    continue
                #A centroidal vector cannot be created for glycine. Throw it out
                if pair[0].residue == "GLY" or pair[1].residue == "GLY":
                    self.message += f"A potential e-g pairing ({pair[0].residue} {pair[0].num} {pair[0].chain}," \
                                    f" {pair[1].residue} {pair[1].num} {pair[1].chain}) contains glycine or the provided" \
                                    f"protein file contains only an alpha carbon. The pair will not be considered\n"
                    continue
                #Correctly order a new pair and add it to a list
                if pair[0].assignment == 'g':
                    vp = self.VectorPair(pair[1], pair[0])
                    newPair = self.AAPair(pair[1], pair[0])
                else:
                    vp = self.VectorPair(pair[0], pair[1])
                    newPair = self.AAPair(pair[0], pair[1])
                #Value setters
                newPair.setVD(vp[0])
                newPair.setEPoint(vp[1])
                newPair.setGPoint(vp[2])
                self.pairs.append(newPair)
        #Clean out bad pairs and get CA, CB, and Cent distances
        self.FilterPairs()
        self.distances =  {"CA":self.GetCAAverage(),
                "CB":self.GetCBAverage(),
                "Centroid":self.GetCentAverage()}

    #Finds and returns the e-g pair with the smallest distance
    def GetPair(self, residue, exclusion):
        minDist = float('inf')
        minAA = None
        for helix in exclusion:
            for aa in [aa for aa in helix if aa.assignment != residue.assignment]:
                if self.GetSimpleDist(residue.centroid, aa.centroid) < minDist:
                    minDist = self.GetSimpleDist(residue.centroid, aa.centroid)
                    minAA = aa
        return [residue, minAA]

    #Gets distance between two objects
    def GetSimpleDist(self, ob1, ob2):
        try:
            return ((ob1[0] - ob2[0]) ** 2 + (ob1[1] - ob2[1]) ** 2 + (ob1[2] - ob2[2]) ** 2) ** (1 / 2)
        except TypeError:
            return float('inf')

    #Get the range of a helix and associated sequence
    def GetRangeSequence(self, helices, chain):
        self.helices = list()
        for helix in helices:
            addHelix = list()
            for i in range(helix[1][0], helix[1][1] + 1):
                changedFlag = False
                try:
                    heptadAssignemnt = helix[3][i - helix[1][0]]
                except IndexError:
                    #Happens if a chain is not continuously listed in its PDB file
                    break
                try:
                    #Add residue to sequence
                    if chain.chains[helix[0]][i].assignment != "none":
                        raise KeyError
                    chain.chains[helix[0]][i].assignment = heptadAssignemnt
                    resToAdd = chain.chains[helix[0]][i]
                    changedFlag = True
                except KeyError:
                    #Happens if a chain is not continuously listed or in some other strange cases.
                    #Iterate through the pdb file until the correct residue is found and label it
                    num = i
                    while num < helix[1][1] + 1:
                        try:
                            num += 1
                            if chain.chains[helix[0]][num].assignment == "none":
                                chain.chains[helix[0]][num].assignment = heptadAssignemnt
                                resToAdd = chain.chains[helix[0]][num]
                                changedFlag = True
                                break
                        except KeyError:
                            continue
                if heptadAssignemnt == 'e' or heptadAssignemnt == 'g' and changedFlag:
                    addHelix.append(resToAdd)
                if resToAdd.num == helix[1][1]:
                    break
            self.helices.append(addHelix)

    #Get vector information for both residues
    #Returns the z-val differernce from the XY orthoginal projection of both vectors,
    #xyz coordinates of the 3D recast for the first residue,
    #and the xyz coordinates of the 3D recast for the second residue
    def VectorPair(self, res1, res2):
        # Slope of line made from CA --> CENTROID
        CA1 = self.GetCA(res1)
        CA2 = self.GetCA(res2)
        # XY orthogonal projection
        m1 = (res1.centroid[1] - CA1.location[1]) / (res1.centroid[0] - CA1.location[0])
        m2 = (res2.centroid[1] - CA2.location[1]) / (res2.centroid[0] - CA2.location[0])
        # Intersection point
        #
        #               y_2 - y_1 + m_1x_1 - m_2x_2
        #   x_int =    ------------------------------
        #                       m_1 - m_2

        x_int = (res2.centroid[1] - res1.centroid[1] + m1 * res1.centroid[0] - m2 * res2.centroid[0]) / (m1 - m2)
        # Cast back into 3D
        #
        #        x - x_1
        #   z = --------- * c + z_1
        #           a

        y1 = (x_int - res1.centroid[0]) / res1.vector[0] * res1.vector[1] + res1.centroid[1]
        y2 = (x_int - res2.centroid[0]) / res2.vector[0] * res2.vector[1] + res2.centroid[1]

        z1 = (x_int - res1.centroid[0]) / res1.vector[0] * res1.vector[2] + res1.centroid[2]
        z2 = (x_int - res2.centroid[0]) / res2.vector[0] * res2.vector[2] + res2.centroid[2]
        return abs(z1 - z2), [x_int, y1, z1], [x_int, y2, z2]

    #Remove unlikely or impossible pairings
    def FilterPairs(self):
        newList = list()
        for val in self.pairs:
            if self.pairs.count(val) > 1 and val not in newList:
                if self.CheckVectorAngle(val) and self.CheckVectorDirection(val):
                    newList.append(val)
        self.pairs = newList

    #Return true if angle <CA1,Cent1,Cent2 angle is acute. Otherwise return false
    def CheckVectorAngle(self, pair):
        ca1 = self.GetCA(pair.e)
        ca2 = self.GetCA(pair.g)
        # <<CA1,CA2,Cent2
        lA1 = self.GetSimpleDist(ca1.location, ca2.location)
        lB1 = self.GetSimpleDist(ca2.location, pair.g.centroid)
        lC1 = self.GetSimpleDist(ca1.location, pair.g.centroid)
        theta1 = math.acos((lA1 ** 2 + lB1 ** 2 - lC1 ** 2) / (2 * lA1 * lB1))
        # <<CA2,CA1,Cent1
        lA2 = self.GetSimpleDist(ca2.location, ca1.location)
        lB2 = self.GetSimpleDist(ca1.location, pair.e.centroid)
        lC2 = self.GetSimpleDist(ca2.location, pair.e.centroid)
        theta2 = math.acos((lA2 ** 2 + lB2 ** 2 - lC2 ** 2) / (2 * lA2 * lB2))
        if math.degrees(theta1) > 90 or math.degrees(theta2) > 90:
            return False
        return True

    #Return true if vectors does not need to be reversed to point towards its xyz intersection point
    def CheckVectorDirection(self, pair):
        eVect = [pair.EPoint[0] - pair.e.centroid[0], pair.EPoint[1] - pair.e.centroid[1],
                 pair.EPoint[2] - pair.e.centroid[2]]
        gVect = [pair.GPoint[0] - pair.g.centroid[0], pair.GPoint[1] - pair.g.centroid[1],
                 pair.GPoint[2] - pair.g.centroid[2]]
        if (np.signbit(eVect).tolist() == (~np.signbit(pair.e.vector)).tolist()) or (np.signbit(gVect).tolist() == (~np.signbit(pair.g.vector)).tolist()):
            return False
        return True

    #Get alpha carbon
    def GetCA(self, res):
        for atom in res.atoms:
            if atom.id == 'CA':
                return atom

    #Get beta carbon
    def GetCB(self, res):
        for atom in res.atoms:
            if atom.id == 'CB':
                return atom

    #Get average CA distance for all pairs in a protein
    def GetCAAverage(self):
        n = 0
        dist = 0
        for pair in self.pairs:
            dist += self.GetSimpleDist(self.GetCA(pair.e).location, self.GetCA(pair.g).location)
            n += 1
        if n != 0:
            return dist / n
        else:
            return " "

    #Get average CB distance for all pairs in a protein
    def GetCBAverage(self):
        n = 0
        dist = 0
        for pair in self.pairs:
            dist += self.GetSimpleDist(self.GetCB(pair.e).location, self.GetCB(pair.g).location)
            n += 1
        if n != 0:
            return dist / n
        else:
            return " "

    #Get average centroid distance for all pairs in a protein. Note that proteins generated from NMR data are often
    #are often missing atoms, so an accurate average centroid distance may not be possible to obtain
    def GetCentAverage(self):
        n = 0
        dist = 0
        for pair in self.pairs:
            dist += self.GetSimpleDist(pair.e.centroid, pair.g.centroid)
            n += 1
        if n != 0:
            return dist / n
        else:
            return " "