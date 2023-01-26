import numpy as np
from numba import njit, jit

PI = np.pi



class Polymer():
    """
    Class Polymer
    ------------------

    

    Attributes:
    -----------
    
    """

    def __init__(self, N, init_config = 'linear'):
        
        # Initialize
        self.N = N
        self.r = np.zeros((self.N,3))
                
        self.dr = np.zeros((N,3))
        self.ds = np.zeros((N))
        self.c  = np.zeros((N))
    
        
        self.init_tree_foil_polymer()
    
    
    
    
    
    def get_geometry(self):
    """
    
    """
        # Find the orientation and length of segments
        for i in range(0,self.N):
            if i==self.N-1:
                self.dr[i] = self.r[0]-self.r[N-1]
                self.ds[i] = np.linalg.norm(self.dr[i,:])
            else:
                self.dr[i] = self.r[i+1]-self.r[i]
                self.ds[i] = np.linalg.norm(self.dr[,:])
                
        for i in range(1, self.N):
            self.c[i] = (dr[i]-dr[i-1]/ds[i]
            
    
    
    def calc_writhe(self):
        """
            Determine the writhe of the polymer following the double
            integral from White's theorem.
        """
        Wr = 0
        
        for i in range(0, self.N+1):    
            ii = i
            ti = self.dr[i,:]
            if i == self.N:
                ii = 0
                ti = dr[self.N-1,:]
                
            for j in range(0, self.N+1):
                jj = j
                tj = self.dr[j,:]
                if j == self.N:
                    jj = 0
                    tj = self.dr[self.N-1,:]
        
                if ii == jj: 
                    continue

                rij = self.r[ii,:] - self.r[jj,:]
                s = np.linalg.norm(rij)
                Wr += np.dot(np.cross(ti,tj), rij)/s**3

        Wr /= (4*PI)
        
        return Wr
    
    
    
    
    # Initialization geometries
    def init_tree_foil_polymer(self):
    """
        Creates a distribution of points arrange to produce 
        a tree-foil knot
    """
    
    t = np.linspace(0,2*PI, self.N)
    for i in range(0,self.N):
        self.r[i,0] = np.sin(t[i]) + 2*np.sin(2*t[i])
        self.r[i,1] = np.cos(t[i]) - 2*np.cos(2*t[i])
        self.r[i,2] = -np.sin(3*t[i])
        
        
    
