import numpy as np

from modules.polymer import Polymer
from modules.montecarlo import MonteCarlo
from modules.exporter import Exporter




if __name__ == "__main__":
    N = 20
    na = Polymer(N)
    mc = MonteCarlo(na)
    """
    for i in range(0,N):
        print(1/na.c[i])
    """
    
    
    mc.rot_section(4, 9, np.pi/2)
    na.get_geometry()
    #print(na.calc_writhe())
    #print(na.unknoted())
    exporter = Exporter(na, "./")
    exporter.save_XYZ_cylinders(0)
    #exporter.save_XYZ(0)
