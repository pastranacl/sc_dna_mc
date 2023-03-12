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
        Checked! Nonetheless, further confirmations necessary
    """
    
    
    def save_XYZ(self, fileid):
        """
            Save file in XYZ for representation with OVITO
        """
        np.savetxt(self.path + str(fileid) + "_na.xyz", 
                    self.polymer.r, 
                    delimiter="\t",
                    header=str(self.N) + "\n",
                    comments='')
    
    
    def save_XYZ_spheres(self, fileid):
        """
            Save extended file in XYZ for representation with OVITO
        """
        
        filexyz = open(self.path + str(fileid) + "_na.xyz","w")
        filexyz.writelines(str(self.N) + "\n")
        filexyz.writelines("Properties=pos:R:3:aspherical_shape:R:3:color:R:3" + "\n")
        
        for i in range(0,self.N):        
            cm = self.polymer.r[i,:] + self.polymer.ds[i]*self.polymer.t[i,:]/2
            
            filexyz.writelines(str(cm[0]) + "\t" + str(cm[1]) + "\t" + str(cm[2]) +  "\t" +                                 # Center of mass of the cylinder
                               str(self.polymer.rpol) + "\t" + str(self.polymer.rpol) + "\t" + str(self.polymer.rpol) + "\t")              # Dimensions
                               
            # Mark in color regions of curvature different than zero                              
            if self.polymer.c0[i]>0: 
                filexyz.writelines("0.7\t0.0\t0.0\n")
            else:
                filexyz.writelines("0.2\t0.2\t0.2\n")

        filexyz.close() 
        
        
        
    def save_XYZ_cylinders(self, fileid):
        """
            Save extended file in XYZ for representation with OVITO
            
            The file is composed by 10 columns. The 3 first are the 
            coordinates of the center of mass of the cylinder. The next 4 
            are a quaternion (3 coordinates for an axis  of rotation + angle). 
            The next 3 are associated to aspherical_shape dictates the radius,
            asphericity (0 for a cylinder), and the length of the cylinder.
            
            See https://www.ovito.org/forum/topic/quaternion-representation-of-aspherical-particels/
            for a description of the quaternion approximation
        """
        
        filexyz = open(self.path + str(fileid) + "_na.xyz","w")

        filexyz.writelines(str(self.N) + "\n")
        filexyz.writelines("Properties=pos:R:3:orientation:R:4:aspherical_shape:R:3:color:R:3" + "\n")
        
        q = np.zeros((4))
        for i in range(0,self.N):
            
            # Quaternion for the rotation with respect to z axis
            q[0:3] = np.cross(np.array([0,0,1]), self.polymer.t[i,:])
            q[3] = 1 + self.polymer.t[i,2] #q[3] = 1 + dot(v1,v2);
            q /= np.linalg.norm(q)
                        
            cm = self.polymer.r[i,:] + self.polymer.ds[i]*self.polymer.t[i,:]/2
            
            filexyz.writelines(str(cm[0]) + "\t" + str(cm[1]) + "\t" + str(cm[2]) +  "\t" +                                 # Center of mass of the cylinder
                               str(q[0]) + "\t" + str(q[1]) + "\t" + str(q[2]) +  "\t" + str(q[3]) + "\t" +                 # Quaternion for orientation
                               str(self.polymer.rpol) + "\t" + str(0) + "\t" + str(self.polymer.ds[i]) + "\t")              # Dimensions
                               
            # Mark in color regions of curvature different than zero                              
            if self.polymer.c0[i]>0: 
                filexyz.writelines("0.7\t0.0\t0.0\n")
            else:
                filexyz.writelines("0.2\t0.2\t0.2\n")

        filexyz.close() 

