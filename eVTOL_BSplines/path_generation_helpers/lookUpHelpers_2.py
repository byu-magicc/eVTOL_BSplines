#this file implements the creation of look up tables for various degree polynomial and length M polynomial
import numpy as np
import os, sys
from pathlib import Path



from eVTOL_BSplines.path_generation_helpers.control_points_matrix import create_W_Matrix
from eVTOL_BSplines.path_generation_helpers.matrix_helpers import B_init_final

#defines the file locations for a bunch of the tables etc.



#brainstorming here: what should the form of the key be? And what should the naming convention be?
#d or degree?
#l or ell?
#d, M, ell
#d, ell, M
#M, ell, d
#M, d, ell
#ell, d, M
#ell, M, d

#let's go with:
#M, d, l



#defines the npz file for the B matrix itself

#B matrix section
BMatrixFileLocation = 'lookUpTables2/B_matrix.npz'
U1MatrixFileLocation = 'lookUpTables2/U1_matrix.npz'
U2MatrixFileLocation = 'lookUpTables2/U2_matrix.npz'
SigmaMatrixFileLocation = 'lookUpTables2/Sigma_matrix.npz'
VtMatrixFileLocation = 'lookUpTables2/Vt_matrix.npz'
PseudoinverseMatrixFileLocation = 'lookUpTables2/Pseudoinverse_matrix.npz'

#S matrix section
SMatrixFileLocation = 'lookUpTables2/S_matrix.npz'



#creates the class that generates the look up tables for the 
class lookUpTablesGenerator:

    #creates the initialization function
    def __init__(self,
                 M_maximum: int = 100, #the maximum number of intervals of interest
                 highestD: int = 5): #the highest degree of which to worry about

        self.M_maximum = M_maximum

        self.highestD = highestD


    

    #creates the function to generate the B tables, and all such variations
    def generateBEndTables(self,
                           BMatrixFileLocation: str = BMatrixFileLocation,
                           U1MatrixFileLocation: str = U1MatrixFileLocation,
                           U2MatrixFileLocation: str = U2MatrixFileLocation,
                           SigmaMatrixFileLocation: str = SigmaMatrixFileLocation,
                           VtMatrixFileLocation: str = VtMatrixFileLocation,
                           PseudoinverseMatrixFileLocation: str = PseudoinverseMatrixFileLocation):
        
        #gets the full file path for the output files of each respective matrix
        temp1 = os.fspath(Path(__file__).parents[0])
        BMatrixFullFilePath = os.path.abspath(os.path.join(temp1, BMatrixFileLocation))
        U1MatrixFullFilePath = os.path.abspath(os.path.join(temp1, U1MatrixFileLocation))
        U2MatrixFullFilePath = os.path.abspath(os.path.join(temp1, U2MatrixFileLocation))
        SigmaMatrixFullFilePath = os.path.abspath(os.path.join(temp1, SigmaMatrixFileLocation))
        VtMatrixFullFilePath = os.path.abspath(os.path.join(temp1, VtMatrixFileLocation))
        PseudoinverseMatrixFullFilePath = os.path.abspath(os.path.join(temp1, PseudoinverseMatrixFileLocation))
    


        #creates the B dictionary 
        B_dictionary = {}
        #creates the U1 dictionary
        U1_dictionary = {}
        #U2 dictionary
        U2_dictionary = {}
        #Sigma dictionary
        Sigma_dictionary = {}
        #Vt dictionary
        Vt_dictionary = {}
        #pseudoinverse dictionary
        pseudoinverseDictionary = {}

        #iterates over the degrees from 1 to highest degree inclusive
        for d in range(1, (self.highestD + 1)):
            
            #iterates M from degree plus one to M maximum
            for M in range((d+1), self.M_maximum):
                
                #gets the C concatenated together
                B_cat = B_init_final(degree=d,
                                     alpha=1.0,
                                     M=M)
                
                #gets the Singular value Decomposition of the B matrix
                U, S, Vt = np.linalg.svd(B_cat)


                #then, we get the partitions of U
                U1 = U[:,:(2*d)]
                U2 = U[:,(2*d):]

                #gets the Sigma as the diagonal array
                Sigma = np.diag(S)

                #section to ge the moore penrose pseudoinverse
                V = np.transpose(Vt)

                U1_T = np.transpose(U1)

                Sigma_inv_items = []
                #gets the reciprocal of Each item in S
                for item in S:
                    temp = 1.0/item
                    #saves this to the items list
                    Sigma_inv_items.append(temp)
                
                #turns the Sigma inverse items into a diagonal matrix
                Sigma_inv = np.diag(np.array(Sigma_inv_items))
                
                #actually gets the peudoinverse
                B_cat_pseudoinverse = V @ Sigma_inv @ U1_T

                #creates the main key
                mainKey = f"M{M}_d{d}"
                
                #adds each component to its respective dictionary
                B_dictionary[mainKey] = B_cat
                
                U1_dictionary[mainKey] = U1

                U2_dictionary[mainKey] = U2

                Sigma_dictionary[mainKey] = Sigma

                Vt_dictionary[mainKey] = Vt

                pseudoinverseDictionary[mainKey] = B_cat_pseudoinverse

        #writes each of those out to their respective dictionaries
        np.savez(BMatrixFullFilePath, **B_dictionary)
        np.savez(U1MatrixFullFilePath, **U1_dictionary)
        np.savez(U2MatrixFullFilePath, **U2_dictionary)
        np.savez(SigmaMatrixFullFilePath, **Sigma_dictionary)
        np.savez(VtMatrixFullFilePath, **Vt_dictionary)
        np.savez(PseudoinverseMatrixFullFilePath, **pseudoinverseDictionary)
        
        
        
        pass
    

    #defines the function to generate the S lookup tables
    def generateSLookupTables(self,
                              SfileLocation: str = SigmaMatrixFileLocation):


        temp1 = os.fspath(Path(__file__).parents[0])
        intputFilePath = os.path.abspath(os.path.join(temp1, SfileLocation))
        
        #creates the S dictionary
        S_dictionary = {}

        #iterates over each degree
        for d in range(self.highestD):

            #instantiates the W matrix lookUpTablesGenerator
            temp_W_gen = create_W_Matrix(d=d)

            #sets the minimum M, which is always greater than the degree by at least one
            M_min = degree + 1
            
            #iterates over the M's
            for M in range(M_min, (self))




#creates the class that reads everything
class lookUpTableReader:

    #defines the init function
    def __init__(self):

        pass
    


    ###########################################################################################
    #B Section of Junk
    #defines the functions to load all of the B matrices
    def loadBLookupTables(self,
                          BMatrixFileLocation: str = BMatrixFileLocation,
                          U1MatrixFileLocation: str = U1MatrixFileLocation,
                          U2MatrixFileLocation: str = U2MatrixFileLocation,
                          SigmaMatrixFileLocation: str = SigmaMatrixFileLocation,
                          VtMatrixFileLocation: str = VtMatrixFileLocation,
                          PseudoinverseMatrixFileLocation: str = PseudoinverseMatrixFileLocation):
        
        #gets the full file path for the output files of each respective matrix
        temp1 = os.fspath(Path(__file__).parents[0])
        BMatrixFullFilePath = os.path.abspath(os.path.join(temp1, BMatrixFileLocation))
        U1MatrixFullFilePath = os.path.abspath(os.path.join(temp1, U1MatrixFileLocation))
        U2MatrixFullFilePath = os.path.abspath(os.path.join(temp1, U2MatrixFileLocation))
        SigmaMatrixFullFilePath = os.path.abspath(os.path.join(temp1, SigmaMatrixFileLocation))
        VtMatrixFullFilePath = os.path.abspath(os.path.join(temp1, VtMatrixFileLocation))
        PseudoinverseMatrixFullFilePath = os.path.abspath(os.path.join(temp1, PseudoinverseMatrixFileLocation))
    

        #loads each of those related submatrices
        B_loaded = np.load(BMatrixFullFilePath)
        U1_loaded = np.load(U1MatrixFullFilePath)
        U2_loaded = np.load(U2MatrixFullFilePath)
        Sigma_loaded = np.load(SigmaMatrixFullFilePath)
        Vt_loaded = np.load(VtMatrixFullFilePath)
        Pseudoinverse_loaded = np.load(PseudoinverseMatrixFullFilePath)

        #and then turns them back into dictionaries to load temporarily
        self.B_data = {B_key: B_loaded[B_key] for B_key in B_loaded.files}       
        self.U1_data = {U1_key: U1_loaded[U1_key] for U1_key in U1_loaded.files}       
        self.U2_data = {U2_key: U2_loaded[U2_key] for U2_key in U2_loaded.files}       
        self.Sigma_data = {Sigma_key: Sigma_loaded[Sigma_key] for Sigma_key in Sigma_loaded.files}       
        self.Vt_data = {Vt_key: Vt_loaded[Vt_key] for Vt_key in Vt_loaded.files}       
        self.Pseudoinverse_data = {Pseudoinverse_key: Pseudoinverse_loaded[Pseudoinverse_key] for Pseudoinverse_key in Pseudoinverse_loaded.files}       
        
        #returns the things we just found
        return self.B_data, self.U1_data, self.U2_data, self.Sigma_data, self.Vt_data, self.Pseudoinverse_data



    #creates the function to return a specific B matrix
    def getIndividualB(self, M: int, d: int):

        #creates the key
        key = f'M{M}_d{d}'

        return (self.B_data)[key]

    #gets individual U1 matrix
    def getIndividualU1(self, M: int, d: int):


        key = f'M{M}_d{d}'

        return (self.U1_data)[key]


    #gets individual U1 matrix
    def getIndividualU2(self, M: int, d: int):


        key = f'M{M}_d{d}'

        return (self.U2_data)[key]


    #gets individual U1 matrix
    def getIndividualSigma(self, M: int, d: int):


        key = f'M{M}_d{d}'

        return (self.Sigma_data)[key]


    #gets individual U1 matrix
    def getIndividualVt(self, M: int, d: int):


        key = f'M{M}_d{d}'

        return (self.Vt_data)[key]


    #gets individual U1 matrix
    def getIndividualPseudoinverse(self, M: int, d: int):


        key = f'M{M}_d{d}'

        return (self.Pseudoinverse_data)[key]

    #defunes the function for 
    
