from modules.polymer import Polymer
from modules.exporter import Exporter




if __name__ == "__main__":
    N = 500
    na = Polymer(N)
    
    """
    for i in range(0,N):
        print(1/na.c[i])
    """
    
    #print(na.calc_writhe())
    print(na.unknoted())
    exporter = Exporter(na, "./")
    exporter.save_XYZ_cylinders(0)
