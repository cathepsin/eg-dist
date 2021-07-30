import HeptadRepeat as HR

class EGDist:

    def GetDistances(self, sequence, chain, centroid, heptad):
        print("Here we go!")
        for helix in heptad.heptadInfo:
            for residue in helix[2]:
                print(HR.OneToThree(residue))


    def CheckSequence(self):
        print("stub")
