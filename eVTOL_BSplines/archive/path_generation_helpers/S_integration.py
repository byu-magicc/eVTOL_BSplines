#this file implements the integration for the multiple sections and obtains
#the S_M^d matrix

import numpy as np
import os, sys
from pathlib import Path
import numpy as np


temp1 = os.fspath(Path(__file__).parents[0])
outputFileName = os.path.abspath(os.path.join(temp1, 'lookUpTables/degree_5_integrations.npz'))



class S_integration_naieve:


    def __init__(self,
                 degree: int, #degree of polynomial basis function in question
                 M: int): #the number of sections of interest


        #loads
        loadedDictionary = np.load(outputFileName)
        allIntegrations = [loadedDictionary[f'degree_{i}'] for i in range(len(loadedDictionary.files))]
        
        integrationList = allIntegrations[degree]
        
        
        length = degree + M

        #creates the matrix to store
        self.S_m_d = np.zeros((length, length))

        #gets the number of sections of basis
        numSectionsBasis = int(degree + 1)
        #creates a nested for loop 
        #iterates over the rows
        for i in range(length):
            #iterates over the columns
            for j in range(length):


                ##initializes the 

                #gets the front position of the first function which is i plus d - 1
                function_1_front_pos = int(i -(degree))
                #gets the front position of the second function, whick is j plus d - 1
                function_2_front_pos = int(j - (degree))

                #gets the rear position of the first function and second function
                function_1_rear_pos = function_1_front_pos - (numSectionsBasis)
                function_2_rear_pos = function_2_front_pos - (numSectionsBasis)


                #iterates over both functions, and creates a list of all the 

                #If our basis function starts at the origin, and the degree is 2 (so 3 subsections for that basis function)
                #then the "locations" for each of our sections will be: {0, 1, 2}. 


                function_1_pos = {}
                function_2_pos = {}


                #creates the dictionary for the basis functions at their present positions
                for k in range(numSectionsBasis):
                    
                    #sets the current position for the current function,
                    #which is incremented one by one
                    currentPosition = function_1_front_pos + k
                    key = k
                    function_1_pos[key] = currentPosition

                #repeats the same for the second function
                for k in range(numSectionsBasis):
                    currentPosition = function_2_front_pos + k
                    key = k
                    function_2_pos[key] = currentPosition


                #gets the values for 1 and 2
                values_1 = set(function_1_pos.values())
                values_2 = set(function_2_pos.values())

                #gets the intersection of the values
                overlapping_values = values_1 & values_2


                #creates the matches
                matches = []


                for value in overlapping_values:

                    keys_1 = [k for k, v in function_1_pos.items() if v==value]
                    keys_2 = [k for k, v in function_2_pos.items() if v==value]


                    for k1 in keys_1:
                        for k2 in keys_2:
                            matches.append((value, k1, k2))

                #so, the keys correspond to the subsection of integration for which to integrate

                #this would be like k1 = 1, k2 = 0, corresponds to i_10
                #then the value for the starting position is the value at the beginning of the thing

                #creates the total sum variabel
                sum = 0.0

                #iterates over the items in the matches
                for startValue, k1, k2 in matches:
                    
                    #gets the value to be added
                    addedSectionIntegration = integrationList[k1,k2]
                    
                    #checks if we are withing the bounds of integration
                    lowerBound = 0
                    upperBound = M
                    if lowerBound <= startValue and startValue < upperBound:
                        

                        #if this is true, we add it to the integration

                        sum = sum + addedSectionIntegration

                #saves the fum to the integration matrix
                self.S_m_d[i,j] = sum
                potato = 0



