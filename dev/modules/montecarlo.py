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
    
    
    # FOR TRIALS
    def mcstep_trial(self, id0, did, alpha):
        
        idx = self.indexes_trial(id0, did)
        self.rot_section(idx, alpha)
        self.polymer_trial.get_geometry()
        print(self.check_intersect(idx))
        self.accept_trial_conf()
    
    
    def test_check_intersect(self):
        res = _check_intersect(self.polymer.r, 
                               self.polymer.N, 
                               self.polymer.t,
                               self.polymer.ds,
                               self.polymer.dpol,
                               np.array([0,1]))
        print(res)
    
    
        
    def mcstep(self) -> None:
        """
            Attemp a monte carlo step by rotating a set of 
            verticies around an axis
        """
        
        # Trial move and update geometry
        alpha = 2*self.max_alpha*np.random.rand() - self.max_alpha
        id0 = np.random.randint(0,self.polymer.N)
        did = np.random.randint(2, self.max_length_rot)
        idx = self.indexes_trial(id0, did)
        
        # Rotate
        self.rot_section(idx, alpha)
        self.polymer_trial.get_geometry()
        
        # Make polymer_trial the new configuration
        if self.mctrial(idx) == True:
            self.accept_trial_conf()
            
            
    def mctrial(self, idx) -> np.bool:
        """
            Evaluates if a trial configuration is accepted.
            It checks for no intersection, lack of knotteness and 
            finally it uses the Energy as a Metropolis criterion 
            to accept or reject the trial configuration
        """
        
        # 1. Check for unknotedness
        """
        if self.polymer_trial.is_knotted() == True:
            return False
        """
        # 2. Check that there is not crossing
       
        if self.check_intersect(idx) == True:
            return False
        
        # 3. Metropolis Criterion (energy)
        self.polymer_trial.get_total_energy()
        E_trial = self.polymer_trial.E_tot
        E_curr  = self.polymer.E_tot
        
        if E_trial<E_curr:
            return True
        else:
            DE = E_trial - E_curr
            if np.random.rand() < np.exp(-DE/self.kBT): 
                return True
            else:
                return False
        

    
    def rot_section(self, idx, alpha) -> None:
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
        id0 = idx[0]
        idf = idx[-1]
       
        r_rot = np.copy(self.polymer.r)
        k = self.polymer.r[idf,:] - self.polymer.r[id0,:]
        k /= np.linalg.norm(k)
        
        for i in idx:
            v = self.polymer.r[i,:] - self.polymer.r[id0,:]
            r_rot[i,:] = v*np.cos(alpha) +  np.cross(k,v)*np.sin(alpha) + k*(1-np.cos(alpha))*np.dot(k,v)
            r_rot[i,:] += self.polymer.r[id0,:]
    
        self.polymer_trial.r = copy.deepcopy(r_rot)
       

  
    def check_intersect(self, idx) -> np.bool:
        """
            Check if there is overlap in the polymer after a Monte Carlo step.
            Wrapped to outside function to jit-it.
            
            Input:
               idx = np.array, vertices involved in the rotation
            
            Output: 
                intersect = bool, Returns true if two cylinders 
                            intersect each otherTrue
        """
        return _check_intersect(self.polymer_trial.r, 
                                self.polymer_trial.N, 
                                self.polymer_trial.t,
                                self.polymer_trial.ds,
                                self.polymer_trial.dpol, 
                                idx)
    
    
    def indexes_trial(self, id0, did) -> np.array:
        """
            Creates an array with the indexes f the vertices 
            involved in the rotation, accounting for boundary
            conditions.
            
            Input: 
                id0 = Initial index
                did = Lenght of the selection
            
            Output: 
                ids = numpy array with the indexes as indicated
        """
        ids = np.linspace(id0, id0 + did, did, endpoint=False, dtype = int)
        for i in range(0, did):
            if ids[i]>=self.polymer.N:
                ids[i] = ids[i] - self.polymer.N
        return ids
    
    
    def accept_trial_conf(self) -> None:
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
        



@njit('boolean(int64[::1], int64)')
def isinarr(inarr, val) -> np.bool:
    """
        Numba version of isin in numpy specific for ints
        
        Input:
            inarr = np.array, input array to check
            val   = int, value to test
            
        Output:
            res = bool, return 1 if val is in the input array inarr
                        and 0 otherwise
    """
    res = False
    for i in range(0, len(inarr)):
        if inarr[i] == val:
            res = True
            break
        
    return res  


@njit('boolean(float64[:,::1], int64, float64[:,::1], float64[::1], float64, int64[::1])')
def _check_intersect(r, N, t, ds, dpol, idx) -> np.bool:   
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
                        
    https://www.geometrictools.com/Documentation/DistanceLine3Line3.pdf
    
    """
    EPS = 1e-12
    
    for i in range(0,N):
    #for i in idx:
        ti = t[i]
        for j in range(0,N):
            if i==j: continue
            if j==i+1: continue
            if i==j+1: continue
            if j==0 and i==N-1: continue
            if i==0 and j==N-1: continue
        
            tj = t[j]
            ddr = r[i,:] - r[j,:]
            
            sq_dotp_titj = np.dot(ti,tj)**2
            
            delta1  =  np.dot(ddr,ti) - np.dot(ddr,tj)*np.dot(ti,tj)
            delta1 /= (sq_dotp_titj - 1 + EPS)
            
            delta2 = np.dot(ddr, tj) - delta1*np.dot(ti,tj)
            
            # Check if the cross between both lines occurs around the two cylinders
            if (delta1>=0 and delta1<=ds[i]) and (delta2>=0 and delta2<=ds[j]):
                r1m = r[i,:] + delta1*ti
                r2m = r[j,:] + delta2*tj
                D = np.linalg.norm(r2m-r1m)
               
                # Determine if the distance is less than the specified radius
                if D <= dpol:
                    return True
        
    return False

  
