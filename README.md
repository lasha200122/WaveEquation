# WaveEquation

```python3
from matplotlib import pyplot as plt
import numpy as np

# speed -> speed
# start -> start interval
# end -> end interval
# f -> function for u(x,0)
# g -> function for u_t(x,0)
# time -> time from 0 to T 
# h -> Delta x
# k -> Delta t
# n -> number of intervals for X
# m -> number of intervals for t
# N -> number of iterations for analysis solution


# function for u(x,0)
def f(x):
    return np.sin(x)

#function for u_t(x,0)
def g(x):
    return 0


class String:
    def __init__(self,speed,start,end,f,g,time,n,m,N):
        self.speed = speed #speed
        self.start = start #start position 
        self.end = end #end position
        self.time = time #time
        self.n = n # number of intervals for x
        self.m = m # number of intervals for time
        self.N = N # number of iterations for analysis solution
        self.l = abs(self.end - self.start) # length of string
        self.f = f # function for u(x,0)
        self.g = g # function for u_t(x,0 )
        
        self.h = abs(end - start) / n # delta x
        self.k = abs(time) / m
        self.r = self.speed * self.k / self.h # Constant that must be less 1
        
        self.domain = np.arange(self.start,self.end, self.h) # domain for x axis interval [a,b]
        self.timeAxis = np.arange(0,self.time,self.k) # list for time interval [0, time]
        
        self.gridForAnalysis = np.zeros((m,n)) # grid fo analysis method solution
        self.gridForNumerical = np.zeros((m,n)) # grid for numerical method solution
        
        #### Printing Part
        print("Speed: " + str(self.speed))
        print("Starting position: " + str(self.start))
        print("Ending position: " + str(self.end))
        print("Length of string: " + str(self.l))
        print("Time: " + str(self.time))
        print("Number of intervals for x: " + str(self.n))
        print("Number of intervals for time: " + str(self.m))
        print("delta x: " + str(self.h))
        print("delta t: " + str(self.k))
        print("Constant: " + str(self.r))
        ####
        
        
        
        if self.r >= 1:
            print()
            print("Can't solve this equation because Constan is: " + str(self.r) + "and it must be less than 1.")
        else:
            self.gridForAnalysis = self.AnalysisSolution() # getting grid for Analysis solution
            self.gridForNumerical = self.NumericalSolution() # getting grid for Numerical solution
            # self.ErrorPerFrame()
            print()
            for i in range(len(self.gridForAnalysis)):
                plt.plot(self.domain,self.gridForAnalysis[i])
                plt.plot(self.domain,self.gridForNumerical[i])
                plt.show()
            
        print("Done.")
        
           
    def SinT(self,N,t):
        return np.sin(N* np.pi * self.speed * t / self.l)
    
    def CosT(self,N,t):
        return np.sin(N* np.pi * self.speed * t / self.l)
    
    def SinX(self,N):
        return np.sin(N* np.pi * self.domain / self.l)
    
    def NumecricalIntegrationForB(self,N):
        constant = 2 / self.l
        sum = 0
        for i in range(len(self.domain)):
            if i == 0 or i == len(self.domain) -1:
                sum += self.f(self.domain[i]) * np.sin(N* np.pi * self.domain[i] / self.l)
            else:
                sum += 2 * self.f(self.domain[i]) * np.sin(N* np.pi * self.domain[i] / self.l)
        return sum
    
    def NumericalIntegrationForBSharp(self,N):
        constant = 2 / (N * self.speed * np.pi)
        sum = 0
        for i in range(len(self.domain)):
            if i == 0 or i == len(self.domain) -1:
                sum += 2 * self.g(self.domain[i]) * np.sin(N* np.pi * self.domain[i] / self.l)
            else:
                sum += self.g(self.domain[i]) * np.sin(N* np.pi * self.domain[i] / self.l)
        return sum
    
    def SummationForAnalysis(self,N,t):
        timePart = self.NumecricalIntegrationForB(N) * self.CosT(N,t) + self.NumericalIntegrationForBSharp(N) * self.SinT(N,t)
        result = timePart * self.SinX(N)
        return result
    
    def AnalysisSolution(self):
        if (self.N < 1):
            print("Analysis solution for this equation can't be solved because N = " + str(self.N))
            return
        grid = np.zeros((self.m,self.n))
        for i, t in enumerate(self.timeAxis):
            u = []
            for n in range(1,self.N+1): 
                if len(u) == 0:
                    u =np.array(self.SummationForAnalysis(n,t))
                else:
                    u += np.array(self.SummationForAnalysis(n,t))
           
            grid[i] = u
        print()
        print("Analysis solution have been found")
        return grid
    
    def NumericalSolution(self):
        grid = np.zeros((self.m,self.n))
        grid[0] = self.f(self.domain)
        for j in range(1,len(grid)- 1):
            for i in range(1,len(grid[j])-1):
                if j==1:
                    grid[j][i] = grid[0][i] + self.k * self.g(self.domain[i])
                else:
                    grid[j][i] = (2 - 2 * np.power(self.r,2)) * grid[j-1][i] + np.power(self.r,2) * (grid[j-1][i+1] + grid[j-1][i-1]) -  grid[j-2][i]
        return grid
    
    def ErrorPerFrame(self):
        for iter in range(len(self.gridForNumerical)):
            difference = self.gridForNumerical[iter] - self.gridForAnalysis[iter]
            error = sum(abs(difference))
            print() 
            print("Frame" + str(iter) + " error : " + str(error))
            print()


String(1,0,np.pi*2,f,g,10,100,200,1)
```



