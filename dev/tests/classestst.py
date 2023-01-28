import copy

class A:
    def __init__(self, N):
        self.N = N
        
        
        
if __name__ == "__main__":
    class1 = A(5)
    class2 = A(2)
    
    print(class1.N)
    print(class2.N)
    print("==========")

    class1 = copy.copy(class2)
    print(class1.N)
    print(class2.N)

    print("==========")
    class1.N = 2
    class2.N = 5
    print(class1.N)
    print(class2.N)
