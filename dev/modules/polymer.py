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
        
        # Parameters input as dictionary (and derived quantities)
        self.kBT = params["kBT"]
        self.p_length = params["pers_length_p"]
        self.C_stiff = params["torsion_stiffness_C"]
        self.rpol = params["polymer_radius"]
        self.R = params["R"]
        self.N = params["N"]
        self.DLk = params["DLk"]
        self.dpol = 2.*self.rpol
      
        
        # Initialize configuration
        self.r = np.zeros((self.N,3))
        self.init_circular_polymer()
        """
        self.init_leminiscata()
        self.init_tree_foil_polymer()
        self.init_figure8_knot()
        """
        
        # Extract energy and geometric parameters
        self.dr = np.zeros((self.N,3))
        self.ds = np.zeros((self.N))
        self.t = np.zeros((self.N,3))
        self.c = np.zeros((self.N))
        self.c0 = np.zeros((self.N))
        self.E_tot = 0
         
        self.get_geometry()
        self.get_total_energy()
    
        # Assing rest curvatures to a region (check if the input is valid)
        if c0 is not None and len(c0)==self.N:
            self.c0 = c0
        elif c0 is not None and len(c0) != self.N:
            print("The input provided for c0 is not valid")
            self.c0 = np.zeros((self.N))
        else:
            self.c0 = np.zeros((self.N))
    
    def __copy__(self):
        polymer_copy = empty_copy(self)
        return polymer_copy
      
    
    def get_total_energy(self):
        Es = self.stretching_energy()
        Eb = self.bending_energy()
        Et = self.torsion_energy()
        self.E_tot = Es + Eb + Et


    def stretching_energy(self):
    #    return np.sum(self.ds)
        return 0

    def bending_energy(self):
        sc = np.sum( (self.c - self.c0)**2)
        return sc*self.p_length


    def torsion_energy(self):
        Wr = self.calc_writhe()
        return self.C_stiff*(self.DLk - Wr)**2
        
    
    
    def is_knotted(self):
        """
            This function uses the Fary-Milnor theorem to check
            if the polymer is knotted. Since the theorem is a 
            sufficient condition, some unknotted configurations
            migth be discarted.
        """
        scfm = np.dot(np.abs(self.c), self.ds)
        if scfm > 4*PI:
            return False
            
        return True
        
    
        
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
            


def empty_copy(obj):
    class Empty(obj.__class__):
        def __init__(self): 
            pass
        
    newcopy = Empty(  )
    newcopy.__class__ = obj.__class__
    return newcopy



@njit('int32(float64[:,:], float64[:,:], int64)')
def _calc_writhe(r, dr, N):
    """
        Calculate the Writhe of the polymer by considering a
        discrete approximation of the double integral.
        
        Input:
            r  =  Coordinates of the polymer nodes
            dr = Direction vector of each edge 
            N  =  Number of nodes
            
        Output:
            Wr = Writhe of the curve
    """
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
