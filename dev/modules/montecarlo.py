import numpy as np
from numba import jit, njit

from polymer import Polymer




Class MonteCarlo():
    """
    Class MonteCarlo
    ------------------

    

    Attributes:
    -----------
    
    """

    def __init__(self, params) -> None:
        """
        """        
        self.kBT = 4.1
        self.r = 1

        
    def mcstep(polymer):
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
        
    
    
    
    
    
    
    
    
    def rot_section(r, id0, idf, alpha):
    """
        Uses Rodrigues' rotation formula to rotate the set
        of verticies between indexes id0 and idf by the 
        angle alpha
    """
        r_rot = np.copy(r)
        k = r[idf,:] - r[id0,:]
        k /= np.linalg.norm(k)
        
        for i in range(id0+1, idf):
            v = r[i,:] - r[id0,:]
            r_rot[i,:] = v*np.cos(alpha) +  np.cross(k,v)*np.sin(alpha) + k*(1-np.cos(alpha))*np.dot(k,v)
            r_rot[i,:] += r[id0,:]

        return r_rot
    
    
    
    
    
    def check_intersect(r, N, eps):
    """
        Check if there is overlap between the curves in
    
        Input: 
                r  = np.array. Coordinates of the array
                N  =
                dcut = Cut off distance between the cylinders
        Output: 
                intersect = bool, Returns true if two cylinders 
                            intersect each other
    """
    
    intersect = False
    for i in range(0,N):
        idp=i+1
        if idp==N:
            idp=0
            
        ti = r[idp,:]-r[i,:]
        ti /= np.linalg.norm(ti)

        for j in range(i+1,N):
        
            jdp=j+1
            if jdp==N:
                jdp=0

            tj = r[jdp,:]-r[j,:]
            tj /= np.linalg.norm(tj)
            dr = r[idp,:] - r[jdp,:]
            
            delta2 = (np.dot(dr,ti)*np.dot(ti,tj) - np.dot(dr,tj)) /        \
                     ( np.dot(ti,tj) - 1)
            
            delta1 = delta2*np.dot(ti,tj) - np.dot(dr, ti)
            
            # Check if the cross between both lines occurs around the two cylinders
            if (delta1>0 and delta1<2*PI/N) and (delta2>0 and delta2<2*PI/N):
                r1m = r[i,:] + delta1*ti
                r2m = r[j,:] + delta2*tj
                d = np.linalg.norm(r2m-r1m)
                
                # Determine if the distance is less than 
                # the specified radius
                if d<eps:
                    intersect = True
                    break
    
