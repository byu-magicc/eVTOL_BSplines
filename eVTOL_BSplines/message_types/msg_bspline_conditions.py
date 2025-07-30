#creates the bspline conditions message

import numpy as np

#This conditions message tells us the condition of a bspline at a particular point in time.
#That is, potentially its position, velocity, acceleration, etc depending on the degree of the spline

#creates the Bspline Conditions class
class MsgBsplineConditions:

    #defines the initialization function
    def __init__(self,
                 degree: int):
        #saves the degree in here
        self.degree = degree

        #creates the conditions list
        self.conditions = []

    #creates the set conditions
    #conditions is in the order of:
    #1. position
    #2. velocity
    #3. acceleration, etc
    def setConditions(self,
                      conditions: list[np.ndarray]):
        
        #gets the length of the conditions,
        #and raises a value error if it's not the same
        conditionsLength = len(conditions)

        if conditionsLength != self.degree:
            raise ValueError(f"Num Conditions {conditionsLength} is not equal to the set degree: {self.degree}")
        
        #goes through and sets each of the conditions
        self.conditions = conditions


