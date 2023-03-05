import numpy as np

from modules.polymer import Polymer
from modules.montecarlo import MonteCarlo
from modules.exporter import Exporter




if __name__ == "__main__":
    
    
    # Define the dictionary for polymer
    N = 10
    polymer_params = {
        "kBT": 4.1,                     # Energy room temperature [pN nm]
        "pers_length_p": 50,            # Persistence length [nm]
        "torsion_stiffness_C": 100,     # Torsion stiffness length, [nm]
        "polymer_radius":(2*np.pi)/N,   # Radius of the polymer (cross-section) [nm]
        "R": 2,                         # Radius of the circular polymer [nm]
        "N": N                         # Number of elements 
    }
    c0 = np.zeros(N)
    c0[0:5] = 0.01
    na = Polymer(polymer_params, c0)
    #print(na.calc_writhe())
    #print(na.unknoted())
    """
    for i in range(0,N):
        print(1/na.c[i])
    """
    
    
    
    # Trial deformation
    mc = MonteCarlo(na)
    mc.mcstep(0, 5, 0.7*np.pi)
    
    
    
    # Save the data
    exporter = Exporter(na, "./")
    exporter.save_XYZ_cylinders(0)
    #exporter.save_XYZ(0)
