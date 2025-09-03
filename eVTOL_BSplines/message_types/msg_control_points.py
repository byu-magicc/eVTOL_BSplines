#this file implements the control points message. 
# Specifically, I have been having troubles
#with communicating control points under the Following circumstances
#So, I have a Series of Safe flight corridors. And in each corridor, I have initial and final conditions
#which define the control points for the bspline in the center. But the end condition for one corridor
#is the start for another. So there are overlapping control points. I need to leverage these 
#control points for my own advantage. But in order to do that, I need a good message class
#that can handle each section, and the total overlapping control point path as well.

import numpy as np


class MSG_Control_Points:

    #defines the initialization function
    #controlPointsArray_list: a list of arrays, where each array is a set of control points
    def __init__(self,
                 degree: int = 3,
                 controlPointsArray_list: list[np.ndarray] = None,
                 controlPointsArray_complete: np.ndarray = None,
                 controlPointsArray_parsed: list[np.ndarray] = None):
        #saves the degree
        self.degree = degree

        self.controlPointsArray_list = controlPointsArray_list
        self.controlPointsArray_complete = controlPointsArray_complete
        self.controlPointsArray_parsed = controlPointsArray_parsed

    
    def setControlPointsArray_list(self,
                                   controlPointsArray_list: list[np.ndarray]):
        
        self.controlPointsArray_list = controlPointsArray_list

    def setControlPointsArray_complete(self,
                                       controlPointsArray_complete: np.ndarray):
        self.controlPointsArray_complete = controlPointsArray_complete

    def setControlPointsArray_parsed(self,
                                     controlPointsArray_parsed: list[np.ndarray]):
        self.controlPointsArray_parsed = controlPointsArray_parsed


    #defines the function to get each of these three
    def getControlPointsArray_list(self):
        return self.controlPointsArray_list
    
    #function to get the controlPoints array complete
    def getControlPointsArray_complete(self):
        return self.controlPointsArray_complete
    
    def getControlPointsArray_parsed(self):
        return self.controlPointsArray_parsed
    
    #function to get the degree
    def getDegree(self):
        return self.degree



    
