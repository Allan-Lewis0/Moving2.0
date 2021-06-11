import math
class Wall():

    def isWall(self,x,y,z):
        # r1=1e-6
        r2=2e-6
        # if (0<=x<=2e-6)& (z**2+y**2>=r1**2):
        #     return 1,r1
        # if (x==2e-6)& (z**2+y**2>r1**2)&(z**2+y**2<=r2**2):
        #     return 2
        if  (z**2+y**2>=r2**2):
            return 1,r2
        else:
            return 0,0

    def Wall_distance(self,x,y,z):
        # if (0<=x<=2e-6)&(z**2+y**2>=1e-6):
        #     return (1e-6)-math.sqrt(y**2+z**2)
        if z**2+y**2<=(2e-6)**2:
            return (2e-6)-math.sqrt(y**2+z**2)
        else:
            return math.sqrt(y**2+z**2)-(2e-6)


