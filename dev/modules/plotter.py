import numpy as np 
import mayavi.mlab as mlab
    
from modules.polymer import Polymer    
    
    
class Plotter3d():
    """
        Class Plotter3d
        ------------------

        

        Attributes:
        -----------
        
    """
    
    def __init__(self, polymer, params) -> None:
        self.polymer = polymer
        self.params = params
        self.N = params["N"]
        self.c0 = params["c0"]
        self.rpol = params["polymer_radius"]
    
    
    def plot(self) -> None:
        
        r_extend = np.zeros((self.N+1, 3))
        r_extend[0:self.N,:] = np.copy(self.polymer.r)
        r_extend[self.N,:] = np.copy(self.polymer.r[0,:])
        
        fig = mlab.figure(bgcolor=(0,0,0))
        #mlab.plot3d(r_extend[:,0], r_extend[:,1], r_extend[:,2],  range(0,self.N+1),  opacity=0.5, tube_radius=2)
        mlab.plot3d(r_extend[:,0], 
                    r_extend[:,1],
                    r_extend[:,2], 
                    tube_radius=self.rpol)
        mlab.show()
