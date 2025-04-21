#this file implements the efficient implementation of the other functions in the legacy matrix helpers file

import numpy as np

#imports the function to evaluate the uniform basis function at a specific set of points
from matrix_helpers import uniform_basis_function_evaluation

#creates the class to sample the basis functions of the 
class basisFunctionSampler:

    #creates the initialization
    def __init__(self,
                 degree: int, #saves the degree of this thing
                 numSamplesPerSection: int): #saves the number of samples per section of the thing
        
        #seaves the degree for this thing
        self.degree = degree

        #saves the number of sections, which is one plus the degree
        self.numSections = self.degree + 1

        #saves the number of samples per section
        self.numSamplesPerSection = numSamplesPerSection

        #gets the Ts
        self.Ts = (1.0/float(numSamplesPerSection))
        

        #creates a list which stores each of the arrays of the sampled basis functions
        self.sampledList = []


        #calls this function automatically to populate the smapled list
        self.sampleBasisFunction()


    #creates the function to sample the basis function
    def sampleBasisFunction(self):




        #iterates over the number of subsections
        for i in range(self.numSections):
            

            #creates the temp list
            tempList = []
            #iterates over the number of samples per section
            for j in range(self.numSamplesPerSection):

                #gets  the current time, which is i plus j times the Ts sample period
                currentTime = float(i) + j*self.Ts

                #gets the sampled value
                sampledValue = uniform_basis_function_evaluation(time=currentTime,
                                                                 degree=self.degree,
                                                                 alpha=1.0)
                

                #saves the sampled value to the temporary list
                tempList.append(sampledValue)

            #converts the temp list into an array
            tempList = np.array(tempList)
            #appends the temp list to the basis samples
            self.sampledList.append(tempList)



#creates a class to write to and read from the 
class sample_and_store_basis_functions:

    #creates the initialization function
    def __init__(self,
                 topDegree: int): #the highest degree of the basis functions from which to sample
        self.topDegree = topDegree


    

