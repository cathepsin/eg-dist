import HeptadRepeat as HR

class EGDist:
    def GetDistances(self, sequence, chain, centroid, heptad):
        print("Here we go!")
        self.GetRangeSequence(heptad.heptadInfo, chain)
        for helix in heptad.heptadInfo:
            for residue in helix[2]:
                print(HR.OneToThree(residue))


    def CheckSequence(self):
        print("stub")

    def AtomDist(self, atom1, atom2):
        return ((atom1.location[0] - atom2.location[0])**2 + (atom1.location[1] - atom2.location[1])**2 + (atom1.location[2] - atom2.location[2])**2)** (1/2)

    def GetRangeSequence(self, helices, chain):
        self.helices = list()
        for helix in helices:
            addHelix = list()
            for i in range(helix[1][0], helix[1][1] + 1):
                addHelix.append(chain.chains[helix[0]][i])
            self.helices.append(addHelix)
