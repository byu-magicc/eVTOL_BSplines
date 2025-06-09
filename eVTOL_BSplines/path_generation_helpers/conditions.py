#this file creates a class which contains the information for conditions for the thing,
#for to make it all work out

import numpy as np

#creates the conditions class, which means that it stores the vector of mechanical derivatives
#for a particular time stamp for the spline. Like position, velocity, acceleration, etc. In a given
#direction for this thing.

#right now, I think that the way this will work is by setting a main initial and Final conditions
#and then the current conditions of the aircraft, which will probably deviate slightly from the desired
#flight path, will become the new initial conditions


class conditions:

    def __init__(self,
                 dimension: int = 2, #sets the dimensionality of the problem to 2
                 numConditions: int = 3): #sets the number of conditions

        #saves the dimension and number of conditions
        self.dimension = dimension
        self.numConditions = numConditions 

        pass




    #defines the function to set the conditions
    #the shape of the conditions array is like:
    #[pos, vel, accel], where pos, vel, accel are column vectors
    def setConditions(self,
                      conditionsList: list[np.ndarray]):
        
        #saves that
        self.conditionsList = conditionsList

        #passes
        pass



        

        
