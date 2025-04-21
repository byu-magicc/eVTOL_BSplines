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
        

        #at the end, we create the 2d array of sampled points
        self.sampledArray = np.ndarray((0,self.numSamplesPerSection))
        for i in range(self.numSections):
            
            #gets the temp vector
            tempVector = (self.sampledList)[i]

            #reshapes the temp vector
            reshapedVector = tempVector.reshape((1,self.numSamplesPerSection))

            #appends to the sampled arrayu
            self.sampledArray = np.concatenate((self.sampledArray, reshapedVector), axis=0)



#creates class to get list of all basis function storages 

class readWriteBasisFunctions:

    #creates the initialization funciton
    def __init__(self,
                 numSamplesPerSection: int, #the number of samples per section to be utilized tiwth this thing.
                 highestDegree: int):#saves the highest degree of bspline from which to sample

        self.highestDegree = highestDegree
        #sets the number of basis types to get
        self.numBasisTypes = self.highestDegree + 1
        #saves the number of samples per section
        self.numSamplesPerSection = numSamplesPerSection

        self.sampledLists = []

        potato = 0


    #defines the function to write the basis functions samples out to an npy or npz file
    def writeToNpz(self, fileName: str):
        #iterates through and gets each of the samples
        for i in range(self.numBasisTypes):
            #gets the basis function sampler for this particular thing
            currentSampler = basisFunctionSampler(degree=i,
                                                  numSamplesPerSection=self.numSamplesPerSection)
            
            #gets the current sampler list
            currentSamplerList = currentSampler.sampledArray

            self.sampledLists.append(currentSamplerList)
        np.savez(fileName, **{f"array_{i}": arr for i,arr in enumerate(self.sampledLists)})

    
    #creates the function to read from Npz
    def readFromNpz(self, filename: str):

        loaded = np.load(filename)
        #obtains the loaded list 
        loaded_list = [loaded[f"array_{i}"] for i in range(self.numBasisTypes)]

        #returns the loaded list
        return loaded_list





#creates the class to create the convolution
class convolveBasisFunction:

    def __init__(self,
                 degree: int): #the degree of the basis polynomial we will be using
        #saves the degree
        self.degree = degree

    #creates the 


