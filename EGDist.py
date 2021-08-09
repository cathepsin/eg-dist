import HeptadRepeat as HR

class EGDist:
    class AAPair:
        def __init__(self, e, g):
            self.e = e
            self.g = g

    def GetDistances(self):
        print("Here we go!")
        self.pairs = list()
        for helix in self.helices:
            exclusionList = [hel for hel in self.helices if hel != helix]
            for residue in helix:
                self.GetPair(residue, exclusionList)


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
            print("")
