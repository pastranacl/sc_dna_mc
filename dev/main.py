import numpy as np

from modules.polymer import Polymer
from modules.montecarlo import MonteCarlo
from modules.exporter import Exporter




if __name__ == "__main__":
    
    
    # Define the dictionary for polymer
    N = 100
    R = 80
    c0 = np.zeros(N)
    c0[0:N//2] = 0.01
    polymer_params = {
        "kBT": 4.1,                             # Energy room temperature [pN nm]
        "pers_length_p": 50,                    # Persistence length [nm]
        "torsion_stiffness_C": 100,             # Torsion stiffness length, [nm]
        "DLk": 1.,                              # Imposed change in the linking number
        "R": R,                                 # Radius of the circular polymer [nm]
        "polymer_radius": (2*R*np.pi)/N,        # Radius of the polymer (cross-section) [nm]
        "N": N,                                 # Number of elements 
        "c0": c0                                # Preferred curvature of each segment
    }
    
    # Define dictionary for Monte Carlo steps
    mc_params = {
        "max_length_rot": 20,                   # Maximum number of segments to rotate`
        "max_alpha": np.pi/2                   # Maximum rotation angle between polymer regions [rads]
    }
    
    
    
    na = Polymer(polymer_params)
    
    """
    print(na.calc_writhe())
    print(na.unknoted())
    for i in range(0,N):
        print(1/na.c[i])
    """
    
    """
        TODO: 
        - Make private attributes explicit with _function
        - Add the types in the returned functions ->
        - Update the torsion stiffness and the stretching stiffness
        - Add the check to evaluate the situation if two edges are parallel (intersection check)
        -
    """   
    # Trial deformation
    mc = MonteCarlo(na, mc_params)
    exporter = Exporter(na, "./output/")
    f = 0
    for i in range(0, 500):
        mc.mcstep()
        if i % 2:
            exporter.save_XYZ_cylinders(f)
            f+=1
            
    
    # Save the data
    #exporter.save_XYZ_cylinders(0)
    #exporter.save_XYZ_spheres(0)
    #exporter.save_XYZ(0)
    
    
    
    
