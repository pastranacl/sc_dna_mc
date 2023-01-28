import numpy as np


class Exporter():
    """
    Class Exporter
    ------------------

    

    Attributes:
    -----------
    
    """
    def __init__(self, polymer, path):
        self.N = polymer.N
        self.polymer = polymer
        self.path = path
        
    """
        TODO: CHECK WHAT HAPPENS IF YOU MODIFY THE POLYMER
        BETTER JUST GIVE THE OBJECT AS ARGUMENT
    """
    def save_XYZ(self, fileid):
        """
            Save file in XYZ for representation with OVITO
        """
        np.savetxt(str(fileid) + "_na.dat", 
                    self.polymer.r, 
                    delimiter="\t",
                    header=str(self.N) + "\n\n")
    
    
    def save_XYZ_cylinders(self, fileid):
        """
        
            Save extended file in XYZ for representation with OVITO
            
            TODO: Add extender description of how is organized the file
            
            See https://www.ovito.org/forum/topic/quaternion-representation-of-aspherical-particels/
            for a description of the quaternion approximation
        """
        
        filexyz = open(self.path + str(fileid) + "_na.xyz","w")

        filexyz.writelines(str(self.N) + "\n")
        filexyz.writelines("Properties=pos:R:3:orientation:R:4:aspherical_shape:R:3" + "\n")
        
        q = np.zeros((4))
        for i in range(0,self.N):
            
            # Quaternion for the rotation with respect to z axis
            q[0:3] = np.cross(np.array([0,0,1]), self.polymer.t[i,:])
            q[3] = 1 + self.polymer.t[i,2]
            q /= np.linalg.norm(q)
                        
            cm = self.polymer.r[i,:] + self.polymer.ds[i]*self.polymer.t[i,:]
            
            filexyz.writelines(str(cm[0]) + "\t" + str(cm[1]) + "\t" + str(cm[2]) +  "\t" +                  # Center of mass of the cylinder
                               str(q[0]) + "\t" + str(q[1]) + "\t" + str(q[2]) +  "\t" + str(q[3]) + "\t" +  # Quaternion for orientation
                               str(0.05) + "\t" + str(0) + "\t" + str(self.polymer.ds[i]) + "\n")            # Dimensions
                            
        
        filexyz.close() 

