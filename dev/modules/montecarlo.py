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
        self.polymer_test = polymer
    
    
    def mcstep(self, id0, idf, alpha):
                
        self.rot_section(id0, idf, alpha)
        self.polymer_test.get_geometry()
                
        condition=True
        if condition:
            self.polymer = copy.copy(self.polymer_test)
            print(self.check_intersect())
            
        
        

    """
    def mcstep():
        E_curr = polymer.E
        
        polymer_test = copy_object
        
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
        
        self.polymer_test.r = copy.copy(r_rot)
        

    #TODO: Note that in the intersection we only need to check the data that has been rotated!
    def check_intersect(self):
        """
            Check if there is overlap between the curves
        
            Input: 
                    r  = np.array. Coordinates of the array
                    N  =
                    dcut = Cut off distance between the cylinders
            Output: 
                    intersect = bool, Returns true if two cylinders 
                                intersect each other
        """

        intersect = False
        for i in range(0, self.polymer_test.N):
            
            idp=i+1
            if idp==self.polymer_test.N: idp=0
                
            ti = self.polymer_test.r[idp,:] - self.polymer_test.r[i,:]
            ti /= np.linalg.norm(ti)

            for j in range(i+1, self.polymer_test.N):
            
                jdp=j+1
                if jdp == self.polymer_test.N: jdp=0


                tj = self.polymer_test.r[jdp,:] - self.polymer_test.r[j,:]
                tj /= np.linalg.norm(tj)
                
                ddr = self.polymer_test.r[i,:] - self.polymer_test.r[j,:]
                
                # TODO: Check how to make that to work
                if np.dot(ti,tj)**2 == 1:
                    if ddr.all() == 0:
                        intersect = True
                    else:
                        intersect = False
                else:
                    delta2 = (np.dot(ddr,tj) - np.dot(ddr,tj)*np.dot(ti,tj)) / \
                             ( np.dot(ti,tj)**2 - 1)
   
                    delta1 = delta2*np.dot(ti,tj) - np.dot(ddr, ti)
                    
                    # Check if the cross between both lines occurs around the two cylinders
                    if (delta1>0 and delta1 < self.polymer_test.ds[i]) and (delta2>0 and delta2 < self.polymer_test.ds[j]):
                        r1m = self.polymer_test.r[i,:] + delta1*ti
                        r2m = self.polymer_test.r[j,:] + delta2*tj
                        d = np.linalg.norm(r2m-r1m)
                        
                        # Determine if the distance is less than 
                        # the specified radius
                        if d<2*self.polymer_test.rpol:
                            intersect = True
                            print(i)
                            print(j)
                            print(d)
                            break

        return intersect
    
    
