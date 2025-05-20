#this file implements the creation of look up tables for various degree polynomial and length M polynomial
import numpy as np
import os, sys
from pathlib import Path



print(Path)

from eVTOL_BSplines.path_generation_helpers.control_points_matrix import create_W_Matrix


class lookUpTableGenerator:

    def __init__(self,
                 M_maximum: int = 100, #the maximum M point. We will start with 100
                 highestDegree: int = 5) -> None:
        
        #saves the highest degree
        self.highestDegree = highestDegree
        
        #saves the M maximum
        self.M_maximum = M_maximum
        
        pass

    #creates the function to generate the S matrices
    def generateSLookupTables(self,
                              fileLocation: str = "lookUpTables/S_Matrices.npz"):
        temp1 = os.fspath(Path(__file__).parents[0])
        intputFilePath = os.path.abspath(os.path.join(temp1, fileLocation))
        #creates the dictionary to store all of the S matrices
        S_dict = {}

        #instantiates the W matrix generator        
        for degree in range(self.highestDegree):
            #instantiates the W matrix generator
            temp_W_gen = create_W_Matrix(d=degree)

            #the minimum boundary for M is d+1, because the requirement
            #is for M to be greater than degree.

            M_min = degree + 1
            #iterates over all of the M's

            for M in range(M_min, (self.M_maximum + 1)):

                #gets the S matrix for the current degree and M_min
                current_C = temp_W_gen.S_k_M(k=degree,
                                            M=M)

                #creates the key identifier for each current C matrix
                key = f"d{degree}_M{M}"
                #stores the current C with the key in the dictionary
                S_dict[key] = current_C

        #saves the complete arrays to an npz file
        np.savez(intputFilePath, **S_dict)


        #passes, and does not return anything
        pass

    #defines the function to generate the list of W matrices, which will
    #become useful to us as time goes on

    def generateWLookupTables(self,
                              fileLocation: str = "lookUpTables/W_d_l_M_Matrices.npz"):

        temp1 = os.fspath(Path(__file__).parents[0])
        intputFilePath = os.path.abspath(os.path.join(temp1, fileLocation))
        
        #crates the dictionary to store all of the W matrices
        W_dict = {}
        
        #iterates through the possible degrees
        #adds one because we need that for the python problem
        for degree in range(self.highestDegree + 1):

            
            #defines the L maximum
            L_max = degree


            #iterates through the various combinations of l
            #l is valid from 1 to L_max inclusive, so I need to do 
            #L_max + 1 because .... Python. I love Python
            for ell in range(1, (L_max + 1)):

                #creates the Minimum Mout
                M_min = degree + 1

                for M in range(M_min, self.M_maximum):
                    

                    #obtains the W temp matrix for this
                    W_temp_creator = create_W_Matrix(d=degree)

                    #obtains the W temp matrix itself
                    W_temp = W_temp_creator.W_d_l_M(d=degree,
                                                    l=ell,
                                                    M=M)


                    #saves the W_temp to the dictionary
                    key = f"degree{degree}_l{ell}_M{M}"

                    #saves the data with the key
                    W_dict[key] = W_temp



        #saves the W_dictionary unpacked to an npz
        np.savez(fileLocation, **W_dict)
        
        #passes and does not return anything
        pass





#creates the main function, which will call the two above functions
def main():

    #instantiates the look up table generator
    generator = lookUpTableGenerator(M_maximum=100,
                                     highestDegree=5)




#calls the main function
if __name__ == "__main__":
    main()




