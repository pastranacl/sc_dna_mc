import numpy as np
from tqdm import tqdm

from modules.polymer import Polymer
from modules.montecarlo import MonteCarlo
from modules.exporter import Exporter


if __name__ == "__main__":
    
    
    # Define the dictionary for polymer
    N = 200
    R = 260
    c0 = np.zeros(N)
    c0[0:N//10] = 1/R
    polymer_params = {
        "kBT": 4.1,                             # Energy room temperature [pN nm]
        "pers_length_p": 50,                    # Persistence length [nm]
        "torsion_stiffness_C": 0,             # Torsion stiffness length, [nm]
        "DLk": 1.,                              # Imposed change in the linking number
        "R": R,                                 # Radius of the circular polymer [nm]
        "polymer_radius": (2*R*np.pi)/N,        # Radius of the polymer (cross-section) [nm]
        "N": N,                                 # Number of elements 
        "c0": c0                                # Preferred curvature of each segment
    }
    
    # Define dictionary for Monte Carlo steps
    mc_params = {
        "max_length_rot": N*0.50,               # Maximum number of segments to rotate`
        "max_alpha": np.pi                      # Maximum rotation angle between polymer regions [rads]
    }
    
    
    """
    print(na.calc_writhe())
    print(na.unknoted())
    for i in range(0,N):
        print(1/na.c[i])
    """
    
    """
        TODO: 
        - Make private attributes explicit with _function
        - Update the torsion stiffness and the stretching stiffness
        - Reorganize the function for checking of intersections
        - Add the check to evaluate the situation if two edges are parallel (intersection check)
    """   
    #  Main routine
    na = Polymer(polymer_params)
    mc = MonteCarlo(na, mc_params)
    exporter = Exporter(na, "./output/")
    f = 0
    for i in tqdm(range(0, 20_000)):
        mc.mcstep()
        if i % 2:
            exporter.save_XYZ_cylinders(f)
            f+=1
            
    
    # Save the data
    #exporter.save_XYZ_cylinders(0)
    #exporter.save_XYZ_spheres(0)
    #exporter.save_XYZ(0)
    
    
    
    
