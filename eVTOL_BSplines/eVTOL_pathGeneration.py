#this file implements the code to implement the eVTOL path generation for the aircraft. 

#In particular, this will perform interpolations after receiving the control points.

#imports the important things
import os, sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from bsplinegenerator.bsplines import BsplineEvaluation

sys.path.insert(0,os.fspath(Path(__file__).parents[1]))

temp1 = os.fspath(Path(__file__).parents[0])
temp2 = os.path.abspath(os.path.join(temp1, 'submodules/path_generator'))
sys.path.insert(0,temp2)
tempPath = sys.path


from eVTOL_BSplines.path_generation_helpers.matrix_helpers import *
from eVTOL_BSplines.path_generation_helpers.conditions_helpers import *

class eVTOL_PathGen:


    #creates the initialization function
    def __init__(self,
                 degree: int, #the degree is the degree of the highest polynomial in the equation here. 
                 M: int, #M is the number of intervals of interest (the time intervals over which the )
                 start_time: float = 0.0, #the start time for this thing
                 dimension: int = 2, #dimension is the number of dimensions we are creating the B Spline within. Usually either 2 or 3
                 alpha: float = 1.0):#alpha is the scaling factor for generating the 
        #stores the degreeinitializeControlPoints
        self.degree = degree
        #stores the M variable
        self.M = M
        #saves the start time
        self.start_time = start_time
        #gets the end time
        self.end_time = self.M*alpha + start_time
        #saves the dimension of the problem
        self.dimension = dimension
        #creates variable to store the time or knot vector
        self.knots = uniform_knot_point_generator(M=M, 
                                                  degree=degree,
                                                  alpha=alpha,
                                                  start_time=start_time)

        #gets the B hat matrix and B hat inverse matrix
        self.B_hat, self.B_hat_inv = B_hat_B_hat_inv(degree=degree,
                                                     alpha=alpha,
                                                     M=M)

        #creates the vector for the control points and the control point set status
        self.controlPoints, self.controlPointSetStatus = initializeControlPoints(dimension=dimension,
                                                                                 degree=degree,
                                                                                 M=M)
    
    #creates function to get all of the control points
    def getControlPoints(self):
        return self.controlPoints

    #creates the function to modify individual control points
    def setControlPoint(self,
                        index: int, #the index of the control points vector to modify
                        controlPoint: np.ndarray): #the control point itself to add
        
        #sets the control point
        (self.controlPoints)[:,index] = controlPoint

        #changes the set status to true
        (self.controlPointSetStatus)[:,index] = True

    #creates the function to set the control points associated with a specific time
    #based on initial conditions
    def setControlPoints_conditions(self,
                                    time: float, #the time at which the initial conditions are being held at 
                                    conditions: conditions_d5): #conditions for the particular point
        
        #checks whether this is an invalid time
        if time < self.start_time or time > self.end_time:
            raise ValueError(f"Invalid time {time}, expected between {self.start_time} and {self.end_time}")
        


        #gets the time index
        timeIndex = int(time)

        #gets the conditions matrix
        conditionsMatrix = conditions.getConditionsMatrix()

        #then we get the control points using equation 23 from Uniform BSplines Paper
        controlPoints = conditionsMatrix @ self.B_hat_inv

        #these index points correspond from (timeIndex) to (timeIndex + d). Though Python's indexing can be weird
        (self.controlPoints)[:,timeIndex:(timeIndex + self.degree)] = controlPoints

        
