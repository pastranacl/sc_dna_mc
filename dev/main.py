import numpy as np
from tqdm import tqdm

from modules.polymer import Polymer
from modules.montecarlo import MonteCarlo
from modules.exporter import Exporter


if __name__ == "__main__":
    
    
    # Define the dictionary for polymer
    N = 50
    R = 130
    c0 = np.zeros(N)
    #c0[0:N//10] = 1/R
    polymer_params = {
        "kBT": 4.1,                             # Energy room temperature [pN nm]
        "pers_length_p": 50,                    # Persistence length [nm]
        "torsion_stiffness_C": 100,             # Torsion stiffness length, [nm]
        "DLk": 10.0,                              # Imposed change in the linking number
        "R": R,                                 # Radius of the circular polymer [nm]
        "polymer_radius": 2,                    # Radius of the polymer (cross-section) [nm]
        "N": N,                                 # Number of elements 
        "c0": c0                                # Preferred curvature of each segment
    }
    
    # Define dictionary for Monte Carlo steps
    mc_params = {
        "max_length_rot": N*0.50,               # Maximum number of segments to rotate`
        "max_alpha": np.pi                     # Maximum rotation angle between polymer regions [rads]
    }
    
    
    """
    print(na.calc_writhe())
    print(na.unknoted())
    for i in range(0,N):
        print(1/na.c[i])
    """
    
    """
        TODO: 
        - Make private attributes explicit with __function
        - Add the check to evaluate the situation if two edges are parallel (intersection check)
        - Intersection for when two edges when they are consecutive
        - 
    """   
    
    
    # Test of intersections
    """
    na = Polymer(polymer_params)
    mc = MonteCarlo(na, mc_params)
    exporter = Exporter(na, "./output/")
    mc.mcstep_trial(0, N//2+1, 0.67*np.pi)
    exporter.save_XYZ_cylinders(0)
    exporter.save_XYZ_spheres(0)
    """
    
    # Test intersections previous execution
    
    r = np.loadtxt("./output/polymer/4999_na.xyz" ,
                   delimiter='\t',
                   usecols=(0,1,2),
                   skiprows=2)
 
    na = Polymer(polymer_params, r)
    mc = MonteCarlo(na, mc_params)
    mc.test_check_intersect()
    """
    from mpl_toolkits import mplot3d
    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.scatter(r[:,0], r[:,1], r[:,2])
    for i in range(1, len(r)):
        ax.plot([r[i,0], r[i-1,0]],
                [r[i,1], r[i-1,1]],
                [r[i,2], r[i-1,2]])
        
    ax.plot([r[0,0], r[len(r)-1,0]],
            [r[0,1], r[len(r)-1,1]],
            [r[0,2], r[len(r)-1,2]])
    plt.show()
    """
    
    #********************************************
    #           M A I N   R O U T I N E         #
    # +++++++++++++++++++++++++++++++++++++++++++
    """
    na = Polymer(polymer_params)
    mc = MonteCarlo(na, mc_params)
    exporter = Exporter(na, "./output/",  False)
    f = 0
    for i in tqdm(range(0, 10_000)):
        mc.mcstep()
        if i % 2:
            exporter.save_XYZ_cylinders(f)
            #exporter.save_XYZ_spheres(f)
            f+=1
    """
    # ---------------------------------------------
    # Save the data
    #exporter.save_XYZ_cylinders(0)
    #exporter.save_XYZ_spheres(0)
    #exporter.save_XYZ(0)
    
    
    
    
