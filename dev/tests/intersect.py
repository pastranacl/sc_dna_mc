import numpy as np


def save_XYZ_cylinders(r0, rf, diam) -> None:
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

        N = len(r0)
        filexyz = open("intersect_test.xyz","w")
        filexyz.writelines(str(N) + "\n")
        filexyz.writelines("Properties=pos:R:3:orientation:R:4:aspherical_shape:R:3:color:R:3\n")
        
        q = np.zeros((4))
        for i in range(0,N):
            
            # Quaternion for the rotation with respect to z axis
            
            t1 = np.array(rf[i,:] - r0[i,:])
            ds = np.linalg.norm(t1)
            t1 /= ds
            cm = r0[i,:] + ds*t1/2
            
            
            q[0:3] = np.cross(np.array([0,0,1]), t1)
            q[3] = 1 + t1[2] #q[3] = 1 + dot(v1,v2);
            q /= np.linalg.norm(q)
                        
            
            
            filexyz.writelines(str(cm[0]) + "\t" + str(cm[1]) + "\t" + str(cm[2]) +  "\t" +      # Center of mass of the cylinder
                               str(q[0]) + "\t" + str(q[1]) + "\t" + str(q[2]) +  "\t" + str(q[3]) + "\t" +                 # Quaternion for orientation
                               str(diam/2) + "\t" + str(0) + "\t" + str(ds) + "\t")              # Dimensions
                               
            # Mark in color regions of curvature different than zero                              
            filexyz.writelines("0.7\t0.0\t0.0\n")


        filexyz.close() 





def check_intersect(r0, rf, diam):
    
    n = len(r0)
    
    intersect = False
    for i in range(0, n):
        
        t1 = rf[i,:] - r0[i,:]
        d1 = np.linalg.norm(t1)
        t1 /= d1
        
        for j in range(i+1,n):
            t2 = rf[j,:] - r0[j,:]
            d2 = np.linalg.norm(t2)
            t2 /= d2
            r12 = r0[i,:]-r0[j,:]
            delta_1 = np.dot(r12,t1)  - np.dot(r12, t2)*np.dot(t1,t2)
            delta_1 /= np.dot(t1,t2)**2 - 1
            
            delta_2 = delta_1*np.dot(t1,t2) + np.dot(r12, t2)
            
            print(delta_1)
            print(delta_2)
            
            if (delta_1 >= 0 and delta_1<=d1) and (delta_2 >= 0 and delta_2<=d2): 
                tr1 = r0[i,:] + t1*delta_1
                tr2 = r0[j,:] + t2*delta_2
                
                distance = np.linalg.norm(tr1 - tr2)
                
                if distance < diam:
                    intersect = True 
                    
    
    return intersect




if __name__ == "__main__":
    
    r0 = np.array( [ [0,0,0],
                     [2,-5.,1.] ])
    
    rf = np.array( [ [0,0,2],
                     [2,5.,1.] ] )
    
    diam = 1.05
    save_XYZ_cylinders(r0, rf, diam)
    
    print(check_intersect(r0, rf, diam))
    
