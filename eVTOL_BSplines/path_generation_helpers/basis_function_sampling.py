#this implements the function to sample the basis functions of the equations.


import numpy as np
from matrix_generators_efficient import basisFunctionSampler, readWriteBasisFunctions
import matplotlib.pyplot as plt
import os, sys
from pathlib import Path


writer = readWriteBasisFunctions(numSamplesPerSection=1000,
                             highestDegree=5)



temp1 = os.fspath(Path(__file__).parents[0])
temp2 = os.path.abspath(os.path.join(temp1, 'lookUpTables/degree_5_basis.npz'))

writer.writeToNpz(temp2)




loadedList = writer.readFromNpz(filename=temp2)

#creates the time array
time = np.linspace(0.0,1.0,1000)


#gets the figure and the axis
figures, axes = plt.subplots(len(loadedList), 1, figsize=(6,8), sharex=True)

#iterates over the list to plot the individual arrays
for i in range(len(loadedList)):
    #gets the current array
    currentBasis = loadedList[i]
    
    #sets the number of lines in the basis function
    numLines = i + 1

    currentAxis = axes[i]
    currentAxis.set_title("Basis function " + str(i))
    for j in range(numLines):
        
        #gets the current array
        currentArray = currentBasis[j,:]
        #plots out from the axis
        currentAxis.plot(time, currentArray, linewidth=2)
    currentAxis.grid(True)

plt.show()

potato = 0