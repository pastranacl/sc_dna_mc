Class MonteCarlo():
    """
    Class Polymer
    ------------------

    

    Attributes:
    -----------
    
    """

    def __init__(self, N) -> None:
        
        # Initialise
        self.N = N
        self.r = init_tree_foil_polymer()
        self.r_current = np.copy(r)
        
    
    
    
    
    def init_tree_foil_polymer():
    """
        
    
    """
    r = np.zeros((N,3))
    t = np.linspace(0,2*PI, N)
    for i in range(0,N):
        r[i,0] = np.sin(t[i]) + 2*np.sin(2*t[i])
        r[i,1] = np.cos(t[i]) - 2*np.cos(2*t[i])
        r[i,2] = -np.sin(3*t[i])
        
    return r 
