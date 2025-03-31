#this file implements the classes to set initial conditions for something

import numpy as np

#creates the class for 5th degree waypoint conditions
class conditions_d5:

    #defines the initialization function
    def __init__(self,
                 pos: np.ndarray, #sets the position for a particular point
                 vel: np.ndarray, #sets the velocity for a particular point
                 accel: np.ndarray,#sets the accel for the point
                 jerk: np.ndarray,#sets the jerk for the point
                 snap: np.ndarray):#sets the snap for the point
        
        #saves them all
        self.pos = pos
        self.vel = vel
        self.accel = accel
        self.jerk = jerk
        self.snap = snap
        #constructs the matrix
        self.conditionsMatrix = np.concatenate((pos,vel,accel,jerk,snap), axis = 1)

    #creates the function for to get the the matrix of the conditions
    def getConditionsMatrix(self):
        #returns the conditions matrix
        return self.conditionsMatrix
