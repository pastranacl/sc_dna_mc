import numpy as np

from modules.polymer import Polymer
from modules.montecarlo import MonteCarlo
from modules.exporter import Exporter




if __name__ == "__main__":
    
    
    # Define the dictionary for polymer
    N = 20
    R = 2
    polymer_params = {
        "kBT": 4.1,                             # Energy room temperature [pN nm]
        "pers_length_p": 50,                    # Persistence length [nm]
        "torsion_stiffness_C": 100,             # Torsion stiffness length, [nm]
        "DLk": 1.,                              # Imposed change in the linking number
        "R": R,                                 # Radius of the circular polymer [nm]
        "polymer_radius": (2*R*np.pi)/N,        # Radius of the polymer (cross-section) [nm]
        "N": N                                  # Number of elements 
    }
    c0 = np.zeros(N)
    c0[0:N//2] = 0.01
    na = Polymer(polymer_params, c0)
    #print(na.calc_writhe())
    #print(na.unknoted())
    """
    for i in range(0,N):
        print(1/na.c[i])
    """
    
    """
        TODO: 
        c0s must be input parameter
        mc_params: dictionary
        make private attributes explicit with _function
        Update the potential
    """   
    # Trial deformation
    mc = MonteCarlo(na)
    for i in range(0, 100):
        mc.mcstep(0, N//2, 0.6*np.pi)
    
    
    # Save the data
    exporter = Exporter(na, "./output/")
    exporter.save_XYZ_cylinders(0)
    #exporter.save_XYZ_spheres(0)
    #exporter.save_XYZ(0)
    
    
    
    
