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

    def __init__(self, polymer, params) -> None:
        """
            params are the range of angles, the maximum extent of rotation, etc
        """        
        self.kBT = polymer.kBT
        self.max_length_rot = params["max_length_rot"]
        self.max_alpha = params["max_alpha"]

        self.polymer = polymer 
        self.polymer_trial = copy.deepcopy(polymer)
    
    """
    # FOR TRIALS
    def mcstep(self, id0, idf, alpha):
                
        self.rot_section(id0, idf, alpha)
        self.polymer_trial.get_geometry()

        condition=True
        if condition:
            self.polymer = copy.copy(self.polymer_trial)
            print(self.check_intersect(id0, idf))
    
    """        
        
        

    def mcstep(self):
        
        
        # Trial move and update geometry
        alpha = self.max_alpha*np.random.rand()
        id0 = np.random.randint(0,self.polymer.N-1)
        idf = id0 + np.random.randint(2, self.max_length_rot)
        
        # ids are periodic and later exchange to have id0 as the smallest
        if idf >= self.polymer.N*np.random.rand(): 
            idf -= self.polymer.N
        
        if idf<id0:
            idf = idf + id0
            id0 = idf - id0
            idf = idf - id0

        # Rotate
        self.rot_section(id0, idf, alpha)
        self.polymer_trial.get_geometry()
        
        # Make polymer_trial the new configuration
        if self.mctrial(id0, idf) == True:
            self.accept_trial_conf()
            print("OK")
           
            
            
            
    def mctrial(self, id0, idf):
        """
            Evaluates if a trial configuration is accepted.
            It checks for no intersection, lack of knotteness and 
            finally it uses the Energy as a Metropolis criterion 
            to accept or reject the trial configuration
        """
        
        # 1. Check unknotedness
        if self.polymer_trial.is_knotted() == False:
            return False
        
        # 2. Check that there is not crossing
        if self.check_intersect(id0, idf) == True:
            return False
        
        # 3. Metropolis Criterion (energy)
        self.polymer_trial.get_total_energy()
        E_trial = self.polymer_trial.E_tot
        E_curr  = self.polymer.E_tot
        
        print(E_curr)
        #print(" - - - - - - -")
        print(E_trial)
        if E_trial<E_curr:
            return True
        else:
            DE = E_trial - E_curr
            print(np.exp( -DE/self.kBT))
            if(np.random.rand() < np.exp( -DE/self.kBT)): 
                return True
            else:
                return False
        

    
    def rot_section(self, id0, idf, alpha):
        """
            Uses Rodrigues' rotation formula to rotate a set
            of vertices 
            
            Input:
                id0 = first index vertex to rotate
                idf = last index vertex to rotate
                alpha = rotation angle
                
            Output:
                polymer_trial.r rotated
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
            Check if there is overlap in the polymer after a Monte Carlo step.
            Wrapped to outside function to jit-it.
            
            Input:
                id0 = first index vertex to rotate
                idf = last index vertex to rotate
            
            Output: 
                intersect = bool, Returns true if two cylinders 
                            intersect each otherTrue
        """
        return _check_intersect(self.polymer_trial.r, 
                                self.polymer_trial.N, 
                                self.polymer_trial.t,
                                self.polymer_trial.ds,
                                self.polymer_trial.dpol, 
                                id0, 
                                idf)
    
    
    def accept_trial_conf(self):
        """
            Copy the main properties of the trial configuration to the
            main configuration
        """
        self.polymer.r = np.copy(self.polymer_trial.r)
        self.polymer.dr = np.copy(self.polymer_trial.dr)
        self.polymer.ds = np.copy(self.polymer_trial.ds) 
        self.polymer.t = np.copy(self.polymer_trial.t)
        self.polymer.c = np.copy(self.polymer_trial.c) 
        self.polymer.E_tot = self.polymer_trial.E_tot
        
        return 0





@njit('boolean(float64[:,::1], int64, float64[:,::1], float64[::1], float64, int64, int64)')
def _check_intersect(r, N, t, ds, dpol, id0, idf):   
    """
        Check if there is overlap in the polymer after a Monte Carlo step
        
        Input: 
                r  = np.array. Coordinates of the array
                N  = Number of points of the polymer
                t  = Tangent vectors per edge
                ds = Length of each edge
                dpol = diameter of the polymer
                id0  = Index of the first vertex that has been modified
                idf  = Index of the last vertex that has been modified
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
