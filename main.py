from prody import *
from pylab import *

ion()

if __name__ == '__main__':
    p38 = parsePDB('1p38')
    p38
    showProtein(p38)
    input("Press enter to continue...")


