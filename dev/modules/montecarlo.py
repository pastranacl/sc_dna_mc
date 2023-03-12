import copy
import numpy as np
from numba import jit, njit

from modules.polymer import Polymer


PI = np.pi


class MonteCarlo():
    """
    Class MonteCarlo
    ------------------

    

    Attributes:
    -----------
    
    
    """

    def __init__(self, polymer) -> None:
        """
            params are the range of angles, the maximum extent of rotation, etc
        """        
        self.kBT = 4.1
        self.r = 1
        
        self.polymer = polymer 
        self.polymer_trial = polymer
    
    
    def mcstep(self, id0, idf, alpha):
                
        self.rot_section(id0, idf, alpha)
        self.polymer_trial.get_geometry()

        condition=True
        if condition:
            self.polymer = copy.copy(self.polymer_trial)
            print(self.check_intersect(id0, idf))
            
        
        

    """
    def mcstep():
        E_curr = polymer.E
        
        polymer_trial = copy_object
        
        # Move
        # Check is not crossing
        p=0
        # Update geometry and everything
        
        # Check unknotedness
        p=0
        
        
        # Calculate energy
        E_move = polymer_mv()
        if E_move<E_curr:
            p=1
        else:
            p = np.exp( -()/self.kBT)
        
        
        return p
        
    """
    
    
    def rot_section(self, id0, idf, alpha):
        """
            Uses Rodrigues' rotation formula to rotate the set
            of verticies between indexes id0 and idf by the 
            angle alpha
        """
        r_rot = np.copy(self.polymer.r)
        k = self.polymer.r[idf,:] - self.polymer.r[id0,:]
        k /= np.linalg.norm(k)
        
        for i in range(id0+1, idf):
            v = self.polymer.r[i,:] - self.polymer.r[id0,:]
            r_rot[i,:] = v*np.cos(alpha) +  np.cross(k,v)*np.sin(alpha) + k*(1-np.cos(alpha))*np.dot(k,v)
            r_rot[i,:] += self.polymer.r[id0,:]
        
        self.polymer_trial.r = copy.copy(r_rot)
       

  
    def check_intersect(self, id0, idf):
        """
            TODO: WRITE DOCUMENTATION
        """
        return _check_intersect(self.polymer_trial.r, 
                                self.polymer_trial.N, 
                                self.polymer_trial.t,
                                self.polymer_trial.ds,
                                self.polymer_trial.dpol, 
                                id0, 
                                idf)
    
    
@njit('boolean(float64[:,::1], int64, float64[:,::1], float64[::1], float64, int64, int64)')
def _check_intersect(r, N, t, ds, dpol, id0, idf):   
    """
        Check if there is overlap between the curves
        
        TODO: WRITE DOCUMENTATION
        Input: 
                r  = np.array. Coordinates of the array
                N  =
                dcut = Cut off distance between the cylinders
        Output: 
                intersect = bool, Returns true if two cylinders 
                            intersect each other
    """

    for i in range(0, N):
        ti = t[i]
        for j in range(id0, idf):

            if i==j: continue
            if j==i+1: continue
            if i==j+1: continue
            if j==0 and i==N-1: continue
            if i==0 and j==N-1: continue
            if i>=id0 and i<idf: continue
            
            tj = t[j]
            ddr = r[i,:] - r[j,:]
            
            sq_dotp_titj = np.dot(ti,tj)**2
            
            # TODO: UPDATE FOR THE CASE OF PARALLEL
            if sq_dotp_titj == 1: continue
         
         
            delta2 = (np.dot(ddr,tj) - np.dot(ddr,ti)*np.dot(ti,tj)) / \
                        (sq_dotp_titj - 1)
            delta1 = np.dot(ddr, ti)-delta2*np.dot(ti,tj)
            
            # Check if the cross between both lines occurs around the two cylinders
            if (delta1>0 and delta1<ds[i]) and (delta2>0 and delta2<ds[j]):
                r1m = r[i,:] + delta1*ti
                r2m = r[j,:] + delta2*tj
                d = np.linalg.norm(r2m-r1m)

                # Determine if the distance is less than the specified radius
                if d < dpol:
                    return True
        
    return False
