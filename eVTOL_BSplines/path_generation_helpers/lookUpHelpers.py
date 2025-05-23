#this file implements the creation of look up tables for various degree polynomial and length M polynomial
import numpy as np
import os, sys
from pathlib import Path



from eVTOL_BSplines.path_generation_helpers.control_points_matrix import create_W_Matrix
from eVTOL_BSplines.path_generation_helpers.matrix_helpers import B_init_final

class lookUpTableGenerator:

    def __init__(self,
                 M_maximum: int = 100, #the maximum M point. We will start with 100
                 highestDegree: int = 5) -> None:
        
        #saves the highest degree
        self.highestDegree = highestDegree
        
        #saves the M maximum
        self.M_maximum = M_maximum
        
        pass



    #defines the function to generate the concatenated B Start and End look up tables
    #and it also saves all of the corresponding singular value Decompositions of the following form:
    #[U1, U2] [Sigma; 0] V_T, and storing as U1, U2, Sigma, V_T. 
    def generateBEndTables(self,
                           fileLocation: str = "lookUpTables/B_Matrices.npz"):

        temp1 = os.fspath(Path(__file__).parents[0])
        intputFilePath = os.path.abspath(os.path.join(temp1, fileLocation))

        #creates the B lookup tables, with the svd data
        B_dict = {}


        
        #iterates over the prospective M's and d's
        for degree in range(1,(self.highestDegree + 1)):

            for M in range((degree+1), self.M_maximum):

                #gets the B_cat
                B_cat = B_init_final(degree=degree,
                                     alpha=1.0,
                                     M=M)
                

                #gets the singular value decomposition of this
                U, S, Vt = np.linalg.svd(B_cat)

                #then, we get the partitions of U
                U1 = U[:,:(2*degree)]
                U2 = U[:,(2*degree):]

                #gets the sigma
                Sigma = np.diag(S)
                                
                #gets the key
                mainKey = f"d{degree}_M{M}"

                #saves the main B matrix
                B_dict[f"{mainKey}_B"] = B_cat


                #saves the B subcomponents
                B_dict[f"{mainKey}_U1"] = U1
                B_dict[f"{mainKey}_U2"] = U2
                B_dict[f"{mainKey}_Sigma"] = Sigma
                B_dict[f"{mainKey}_Vt"] = Vt



        #saves the file to that file path
        np.savez(intputFilePath, **B_dict)

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
        np.savez(intputFilePath, **W_dict)
        
        #passes and does not return anything
        pass




#creates the class to read from and process the look up table files
class lookUpTableReader:

    #creates the init function
    def __init__(self):
        pass

    #creates the function to read the S Lookup tables
    def readSLookupTables(self,
                          fileLocation: str = "lookUpTables/S_Matrices.npz"):
    
        temp1 = os.fspath(Path(__file__).parents[0])
        intputFilePath = os.path.abspath(os.path.join(temp1, fileLocation))

        #gets the loaded S matrix
        S_loaded = np.load(intputFilePath)


        #creates a dictionary from the npz file
        S_data = {key: S_loaded[key] for key in S_loaded.files}

        #saves the S_data to self
        self.S_data = S_data
        #returns the data
        return S_data
    
    #creates the function to get the individual S
    def getIndividualS(self, d: int, M: int):

        key = f"d{d}_M{M}"

        #returns the S for that key
        return (self.S_data)[key]

    #creates the function to read the B Lookup tables
    def readBLookupTable(self,
                         fileLocation: str = "lookUpTables/B_Matrices.npz"):
        
        temp1 = os.fspath(Path(__file__).parents[0])
        inputFilePath = os.path.abspath(os.path.join(temp1, fileLocation))

        #gets the loaded B matrix
        B_loaded = np.load(inputFilePath)

        #creates a dictionary from the npz file
        B_data = {key: B_loaded[key] for key in B_loaded.files}

        #saves the B_data
        self.B_data = B_data
        
        #returns the data
        return B_data
    
    #creates a function to read in the B lookup table dictionary, and get the information
    #corresponding to a particular d and M
    def getIndividualB(self, d: int, M: int):

        #creates the key
        mainKey = f'd{d}_M{M}'

        #gets the B_temp
        B = (self.B_data)[f'{mainKey}_B']
        U1 = (self.B_data)[f'{mainKey}_U1']
        U2 = (self.B_data)[f'{mainKey}_U2']
        Sigma = (self.B_data)[f'{mainKey}_Sigma']
        Vt = (self.B_data)[f'{mainKey}_Vt']

        #returns the components
        return B, U1, U2, Sigma, Vt



    #creates the function to read the W Lookup tables
    def readWLookupTable(self,
                         fileLocation: str = "lookUpTables/W_d_l_M_Matrices.npz"):
        
        temp1 = os.fspath(Path(__file__).parents[0])
        inputFilePath = os.path.abspath(os.path.join(temp1, fileLocation))

        #gets the loaded B matrix
        W_loaded = np.load(inputFilePath)

        #creates a dictionary from the npz file
        W_data = {key: W_loaded[key] for key in W_loaded.files}
        
        #saves the W data
        self.W_data = W_data

        #returns the data
        return W_data
    
    def getIndividualW(self, d: int, ell: int, M: int):

        key = f"degree{d}_l{ell}_M{M}"

        W = (self.W_data)[key]

        return W