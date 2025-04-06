#%%
import os, sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from bsplinegenerator.bsplines import BsplineEvaluation
from IPython.display import display

from IPython.display import display

import scipy.integrate as integrate

sys.path.insert(0,os.fspath(Path(__file__).parents[1]))
tempPath = sys.path


from eVTOL_BSplines.path_generation_helpers.conditions_helpers import conditions, conditionsList

from eVTOL_BSplines.path_generation_helpers.matrix_helpers import b_d_M_t_vector



#'''
#defines the function to get a value the matrix at a given time
def integrationMatrix(t: float, #time (independent variable to integrate over)
                      degree: int, #degree of the bspline
                      alpha: float, #scaling factor
                      l: int, #degree of derivative
                      M: int, #number of integrals of interes
                      row: int, #the row number
                      col: int): #the column number 
    
    #gets the matrix degree
    vectorDegree = degree - l
    vector = b_d_M_t_vector(time=t,degree=vectorDegree,alpha=alpha,M=M)
    #gets the matrix
    A_output = vector @ np.transpose(vector)
    
    #gets the value of the matrix at a particular row and column
    outputValue = A_output[row,col]
    #returns the Matrix
    return outputValue



M = 15
degree = 5
l = 4

#gets the length of the matrix
matrixLength = M + degree - l


#

IntegratedMatrix = np.zeros((matrixLength, matrixLength))
#iterates down through the rows
for i in range(matrixLength):
    #iterates across the row
    for j in range(matrixLength):
        #integrates the components
        tempValue = integrate.quad(func=integrationMatrix, a=0, b=M, args=(degree,1.0,l,M,i,j))
        IntegratedMatrix[i,j] = tempValue[0]


#outputs the matrix
display(IntegratedMatrix)
#'''

'''
def func(x: float, row: int, col: int):

    array = np.array([[1,x],
                      [x**2,x**3]])
    return array[row, col]

outputMatrix = np.zeros((2,2))
for i in range(2):
    for j in range(2):
        tempresult = integrate.quad(func=func, a=0, b=1, args=(i,j))
        outputMatrix[i,j] = tempresult[0]

        banana = 0


display(outputMatrix)
#'''





# %%
