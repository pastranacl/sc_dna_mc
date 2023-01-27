from modules.polymer import Polymer






if __name__ == "__main__":
    N = 500
    na = Polymer(N)
    
    """
    for i in range(0,N):
        print(1/na.c[i])
    """
    
    print(na.calc_writhe())
