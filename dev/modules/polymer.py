import numpy as np
from numba import jit, njit
from numba.experimental import jitclass


PI = np.pi

"""
spec = [
    ('value', int32),               # a simple scalar field
    ('array', float32[:]),          # an array field
]
@jitclass(spec)
"""

class Polymer():
    """
    Class Polymer
    ------------------

    

    Attributes:
    -----------
    
    """

    def __init__(self, params, c0=None):
        
        # Parameters input as dictiontary
        self.kBT = params["kBT"]
        self.p_length = params["pers_length_p"]
        self.C_stiff = params["torsion_stiffness_C"]
        self.rpol = params["polymer_radius"]
        self.R = params["R"]
        self.N = params["N"]
        


        # Initialize
        self.r = np.zeros((self.N,3))
        self.init_circular_polymer()
        #self.init_leminiscata()
        #self.init_tree_foil_polymer()
        #self.init_figure8_knot()
        
        # Extract geometric parameters
        self.dr = np.zeros((self.N,3))
        self.ds = np.zeros((self.N))
        self.t = np.zeros((self.N,3))
        self.c = np.zeros((self.N))
        # self.c0 = c0 if c0 is not None and len(c0)==self.N else np.zeros((self.N))
        if c0 is not None and len(c0)==self.N:
            self.c0 = c0
        elif c0 is not None and len(c0) != self.N:
            print("c0 provided is not valid")
            self.c0 = np.zeros((self.N))
        else:
            self.c0 = np.zeros((self.N))

        self.get_geometry()
    
    
    
    def total_energy(self):
        get_geometry()
        Es = stretching_energy()
        Eb = bending_energy()
        Et = torsion_energy()
        return Es + Eb + Et


    def stretching_energy(self):
        return np.sum(self.ds)
    

    def bending_energy(self):
        sc = np.sum(self.c - self.c0)
        return sc*BF


    def torsion_energy(self):
        Wr = calc_writhe()
        return (DLk - Wr)**2
        
    
    
    def unknoted(self):
        """
            This function uses the Fary-Milnor theorem to check
            if the polymer us unknotted. Since the theorem is a 
            sufficient condition, some unknotted configurations
            migth be discarted.
        """
        uknotted = True
        scfm = np.dot(np.abs(self.c), self.ds)
        if scfm > 4.*PI:
            uknotted = False
        return uknotted
        
    
        
    def get_geometry(self):
        """
            Get the geometry parameters: distance between cylinders
            tangent vectors and curvatures
        """
        
        # i. orientation and length of segments
        for i in range(0,self.N):
            if i==self.N-1:
                self.dr[i,:] = self.r[0,:]-self.r[self.N-1,:]
            else:
                self.dr[i,:] = self.r[i+1,:]-self.r[i,:]
                
            self.ds[i] = np.linalg.norm(self.dr[i,:])
            self.t[i] = self.dr[i,:]/self.ds[i]
            
        # ii. local curvatures
        for i in range(0, self.N):
            if i==0:
                self.c[0] = np.linalg.norm((self.t[0,:] - self.t[self.N-1,:]) / self.ds[0])
            else:
                self.c[i] = np.linalg.norm((self.t[i,:] - self.t[i-1,:]) / self.ds[i])
        


    def calc_writhe(self):
        """
            Determine the writhe of the polymer following the double
            integral from White's theorem.
            This part is computationally expensive since it grows as 
            O(N^2). Thus the function is jited and as such we add a 
            wrapper and put it out of the class.
        """
        Wr = _calc_writhe(self.r, self.dr, self.N)
        return Wr
    
    
    ##########################################################
    #                Initialization geometries               #
    ##########################################################
    
    def init_circular_polymer(self):
        """
            Creates a distribution of points arrange to produce 
            a circular ring of radius 1
        """

        dpolar = 2.0*np.pi/self.N
        for i in range(0, self.N):
            self.r[i,0] = self.R*np.cos(dpolar*i)
            self.r[i,1] = self.R*np.sin(dpolar*i)
            self.r[i,2] = 0


    def init_tree_foil_polymer(self):
        """
            Creates a distribution of points arranged to produce 
            a tree-foil knot (Wr = -3)
        """
        p = np.linspace(0,2*PI, self.N)
        for i in range(0,self.N):
            self.r[i,0] = np.sin(p[i]) + 2*np.sin(2*p[i])
            self.r[i,1] = np.cos(p[i]) - 2*np.cos(2*p[i])
            self.r[i,2] = -np.sin(3*p[i])

        
    def init_figure8_knot(self):
        """
            Creates a distribution of points arranged to produce 
            the Figure 8 knot (bretzel-like) Wr = 0
        """
        p = np.linspace(0,2*PI,self.N)
        for i in range(0, self.N):
            self.r[i,0] = (2 + np.cos(2*p[i]))*np.cos(3*p[i])
            self.r[i,1] = (2 + np.cos(2*p[i]))*np.sin(3*p[i])
            self.r[i,2] = np.sin(4*p[i])

            
    def init_leminiscata(self):
        """
            Creates a distribution of points arranged to produce 
            a leminiscata shape
        """
        dpolar = 2.0*np.pi/self.N
        a=5
        for i in range(0, self.N):
            self.r[i,0] = a*np.cos(dpolar*i)/(1 + np.sin(dpolar*i)**2 )
            self.r[i,1] = a*np.cos(dpolar*i)*np.sin(dpolar*i)/(1 + np.sin(dpolar*i)**2 )
            self.r[i,2] = 0
            






@njit('int32(float32[:,:], float32[:,:], int32)')
def _calc_writhe(r, dr, N):
    Wr = 0
    rij = np.zeros(3)
    for i in range(0,N):    
        ii = i
        ti = dr[i,:]
        if i == N:
            ii = 0
            ti = dr[N-1,:]
            
        for j in range(0,N):
            jj = j
            tj = dr[j,:]
            if j == N:
                jj = 0
                tj = dr[N-1,:]
    
            if ii == jj: 
                continue

            rij = r[ii,:] - r[jj,:]
            s = np.linalg.norm(rij)
            Wr += np.dot(np.cross(ti,tj), rij)/s**3

    Wr /= (4*PI)
    
    return Wr
