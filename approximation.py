import scipy.integrate as integrate
import numpy as np

def approx(arbitraryFunction, xmin, xmax, nWells, Eo, sigma):
    nWells += 1
    xSpace = np.linspace(xmin, xmax, 1000, retstep=True)    #Set of points whithin which the function should be welled, with a set precision
    step = xSpace[1]                                        #linspace gives the step thanks to the retstep argument
    functionMinimum = min(arbitraryFunction(xSpace[0], Eo, sigma))     #The minimum of the function within the desired space
    
    approx = []                                             #Empty list that will be filled with the calculations

    #Offsets the entire arbitrary function to avoid negative integrals (so the function is always positive)
    def f(x):                                               
         return arbitraryFunction(x, Eo, sigma) + np.abs(functionMinimum)
    
    #Dividing the total area calculated with an integral by the number of wells ( that should have an equal area )
    area = np.abs(integrate.quad(f, xmin, xmax))
    wellArea = area[0]/nWells
    
    #Setting the well interval
    tempMin = xmin
    i = tempMin

    #Tracking the current well and point (debug purposes)
    wellNb = 0 
    pointCount = 0
    
    #Giving the approximation an initial value at the leftmost point
    approx.append(f(xmin))
    
    totalPoint = 1
    
    #Well calculation
    while i < xmax:

        localArea = integrate.quad(f, tempMin, i)
        #print("i =", i, "tempMin", tempMin, "xmax", xmax, "localArea", localArea, "wellArea", wellArea, "wellNo", wellNb, "point", pointCount)
        
        if(np.abs(localArea[0]) >= np.abs(wellArea)):
            #print("localArea = ", localArea[0],">= wellArea =",wellArea, "for x range:(", xmin, ";", i, ")",i -xmin)
            j = tempMin
            localInterval = np.arange(j,i,step)
            localMin = min(f(localInterval))
            
            while j < i:
                    totalPoint += 1
                    #print("append " ,i, "f:", localMin, "POINTCOUNT" ,totalPoint)
                    approx.append(localMin)
                    j = j+step
            tempMin = i
            wellNb += 1
            
            i = i+step
    
        #End of interval case
        elif(np.abs(i - xmax) <= step):
                
                j = tempMin + step
                
                while j < xmax:
                    totalPoint +=1
                    #print("append (ELIF)", j, "f:", f(tempMin), "POINTCOUNT" ,totalPoint)
                    approx.append(f(tempMin))
                    j = j+step
                
                break
        else:
                i = i + step
        pointCount += 1

    if(totalPoint < 1000):
          approx.append(f(tempMin))
          print("append needed :c")
    #reverse the initial offset
    return approx - np.abs(functionMinimum)


