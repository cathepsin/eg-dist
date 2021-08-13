import HeptadRepeat as HR
import operator

class EGDist:
    class AAPair:
        def __init__(self, e, g):
            self.e = e
            self.g = g

        def setVD(self, val):
            self.zdist = val

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

        def GetDist(self, a1, a2):
            return ((a1[0] - a2[0]) ** 2 + (a1[1] - a2[1]) ** 2 + (a1[2] - a2[2]) ** 2) ** (1 / 2)

        def CADist(self):
            #Get CAs
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

        def CBDist(self):
            #Get CBs
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

        def CentDist(self):
            return self.GetDist(self.e.centroid, self.g.centroid)



    def GetDistances(self):
        self.pairs = list()
        for helix in self.helices:
            exclusionList = [hel for hel in self.helices if hel != helix]
            for residue in helix:
                pair = self.GetPair(residue, exclusionList)
                if pair[0].assignment == 'g':
                    print(self.pairs.count(self.AAPair(pair[1],pair[0])))
                    self.pairs.append(self.AAPair(pair[1],pair[0]))
                else:
                    print(self.pairs.count(self.AAPair(pair[0], pair[1])))
                    self.pairs.append(self.AAPair(pair[0], pair[1]))
                self.pairs[len(self.pairs) - 1].setVD(self.VectorPair(pair[0], pair[1]))
        self.FilterPairs()



    def GetPair(self, residue, exclusion):
        minDist = float('inf')
        minAA = None
        for helix in exclusion:
            for aa in [aa for aa in helix if aa.assignment != residue.assignment ]:
                if self.CentroidDist(residue.centroid, aa.centroid) < minDist:
                    minDist = self.CentroidDist(residue.centroid, aa.centroid)
                    minAA = aa
        return [residue, minAA]

    def CentroidDist(self, cen1, cen2):
        return ((cen1[0] - cen2[0])**2 + (cen1[1] - cen2[1])**2 + (cen1[2] - cen2[2])**2)** (1/2)

    def GetRangeSequence(self, helices, chain):
        self.helices = list()
        for helix in helices:
            addHelix = list()
            for i in range(helix[1][0], helix[1][1] + 1):
                heptadAssignemnt = helix[3][i - helix[1][0]]
                chain.chains[helix[0]][i].assignment = heptadAssignemnt
                if heptadAssignemnt == 'e' or heptadAssignemnt == 'g':
                    addHelix.append(chain.chains[helix[0]][i])
            self.helices.append(addHelix)

    def VectorPair(self, res1, res2):
        # Slope of line made from CA --> CENTROID
        CA1 = None
        for atom in res1.atoms:
            if atom.id == 'CA':
                CA1 = atom
                break
        CA2 = None
        for atom in res2.atoms:
            if atom.id == 'CA':
                CA2 = atom
                break
        #XY orthogonal projection
        m1 = (res1.centroid[1] - CA1.location[1]) / (res1.centroid[0] - CA1.location[0])
        m2 = (res2.centroid[1] - CA2.location[1]) / (res2.centroid[0] - CA2.location[0])
        # Intersection point
        #
        #               y_2 - y_1 + m_1x_1 - m_2x_2
        #   x_int =  --------------------------
        #                       m_1 - m_2

        x_int = (res2.centroid[1] - res1.centroid[1] + m1 * res1.centroid[0] - m2 * res2.centroid[0]) / (m1 - m2)
        print("X intercept of ", res1, " and ", res2, ":", x_int)
        #Cast back into 3D
        #
        #        x - x_1
        #   z = --------- * c + z_1
        #           a
        z1 = (x_int - res1.centroid[0])/res1.vector[0] * res1.vector[2] + res1.centroid[2]
        z2 = (x_int - res2.centroid[0]) / res2.vector[0] * res2.vector[2] + res2.centroid[2]
        return abs(z1 - z2)

    def FilterPairs(self):
        newList = list()
        for val in self.pairs:
            if self.pairs.count(val) > 1 and val not in newList:
                newList.append(val)
        self.pairs = newList
        print("clean")

    def CheckVectorAngle(self, CA_interest, CA_other, x_int):
        #TODO Check the angle <CA1,CA2,x_int. If obtuse, throw out the pairing
        print("stub")

    def CheckVectorDirection(self):
        #TODO Determine if a vector must be multiplied by a negative scaler to reach the x_int. If so, throw out the pairing
        print("stub")