#this file implements the classes to set initial conditions for something

import numpy as np

#creates the class for 5th degree waypoint conditions
class conditions:

    #defines the initialization function
    def __init__(self,
                 numDerivatives: int, #sets the degree of the derivative which to use
                 dimension: int, #sets the dimension of the spline, usually either 2d or 3d
                 time: float): #sets the time of this condition


        #saves the time
        self.time = time
        #saves the number of derivatives
        self.numDerivatives = numDerivatives
        #saves the dimension of the problem
        self.dimension = dimension

        #saves them all
        self.pos = np.zeros((dimension,1))
        self.vel = np.zeros((dimension,1))
        self.accel = np.zeros((dimension,1))
        self.jerk = np.zeros((dimension,1))
        self.snap = np.zeros((dimension,1))
        #constructs the matrix
        self.conditionsMatrix = np.concatenate((self.pos,self.vel,self.accel,self.jerk,self.snap), axis = 1)

    #creates the function to set the position
    def setPosition(self, pos: np.ndarray):
        self.pos = pos
        self.updateConditionsMatrix()

    #creates the function to set the velocity
    def setVelocity(self, vel: np.ndarray):
        self.vel = vel
        self.updateConditionsMatrix()

    #sets acceleration
    def setAccel(self, accel: np.ndarray):
        self.accel = accel
        self.updateConditionsMatrix()

    #sets the jerk
    def setJerk(self, jerk: np.ndarray):
        self.jerk = jerk
        self.updateConditionsMatrix()

    #sets the snap
    def setSnap(self, snap: np.ndarray):
        self.snap = snap
        self.updateConditionsMatrix()


    #creates the update matrix function
    def updateConditionsMatrix(self):
        self.conditionsMatrix = np.concatenate((self.pos,self.vel,self.accel,self.jerk,self.snap), axis = 1)

    #creates the function for to get the the matrix of the conditions
    def getConditionsMatrix(self):
        self.updateConditionsMatrix()
        #gets the section of the conditions matrix to retrieve
        shortenedMatrix = ((self.conditionsMatrix)[:,:(self.numDerivatives+1)]).reshape((self.dimension, (self.numDerivatives+1)))
        #returns the conditions matrix
        return shortenedMatrix
    



#creates the class that contains a list of conditions from the above conditions class
class conditionsList:

    #creates the initialization function
    def __init__(self):

        #creates the array to store the list of conditions
        self.allConditions = []

    #creates a function to add a conditions object
    def addCondition(self, new_condition: conditions):

        #adds the new condition to the array
        (self.allConditions).append(new_condition)

    #creates a function to get the array of all the conditions
    def getAllConditions(self):
        #returns it
        return self.allConditions
    
    #creates a function to get a specific conditions
    def getCondition(self, index: int):

        return (self.allConditions)[index]

    #creates a concatenated matrix of conditions
    def getAllConditionsMatrix(self):
        
        #gets the dimension of the problem
        dimension = self.allConditions[0].dimension

        allConditionsMatrix = np.ndarray((dimension,0))
        #iterates through them all
        for currentCondition in self.allConditions:
            #gets the temp matrix
            tempMatrix = currentCondition.getConditionsMatrix()
            #concatenates it
            allConditionsMatrix = np.concatenate((allConditionsMatrix, tempMatrix), axis=1)

        #returns the all conditions matrix
        return allConditionsMatrix

    


