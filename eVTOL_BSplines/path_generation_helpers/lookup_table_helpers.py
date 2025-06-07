#this file implements the lookup table generators and readers classes
import numpy as np
import os, sys
from pathlib import Path


from eVTOL_BSplines.path_generation_helpers.W_matrix_generator import create_W_Matrix
from eVTOL_BSplines.path_generation_helpers.general_matrix_helpers import B_init_final
from eVTOL_BSplines.path_generation_helpers.basis_function_helpers import readWriteBasisFunctions, integrateBasisFunctionsContinuous


#sets some parameters
#the sample rate per section
numSamplesPerSection = 1000
#the highest degree. We will sample all of the degrees below this
#so, we will choose degree 5, so it will sample degrees 0, 1, 2, 3, 4, and 5 inclusive
highestDegree = 5

#sets the number of intervals of interest
numIntervalsOfInterest = 100


#gets the currentfile path
current_file_path = Path(__file__).resolve()
#gets the project root
innerPackageDirectory = current_file_path.parents[1]
#gets the lookuptable directory
lookupDirectory = innerPackageDirectory / 'lookUpTables'

#defines the npz file for the B matrix itself
#basis function section
#remember, that this is up to degree 5 polynomaials and so forth

#the location of the sampled basis function lookup table
BasisFunctionSampledFileName = 'BasisFunctionSampled.npz'
#the location of the integrated basis function lookup table
BasisFunctionIntegratedFileName = 'BasisFunctionIntegrated.npz'


#B matrix section
BMatrixFileName = 'B_matrix.npz'
U1MatrixFileName = 'U1_matrix.npz'
U2MatrixFileName = 'U2_matrix.npz'
SigmaMatrixFileName = 'Sigma_matrix.npz'
VtMatrixFileName = 'Vt_matrix.npz'
PseudoinverseMatrixFileName = 'Pseudoinverse_matrix.npz'

#S matrix section
SMatrixFileName = 'S_matrix.npz'

#W matrix section
WMatrixFileName = 'W_matrix.npz'

#Y and Z file locations
YMatrixFileName = "Y_Matrices.npz"
ZMatrixFileName = "Z_Matrices.npz"


#gets the full directories for the matrices and everything, using the lookup Directory

#Basis function section
BasisFunctionSampledDirectory = lookupDirectory / BasisFunctionSampledFileName
BasisFunctionIntegratedDirectory = lookupDirectory / BasisFunctionIntegratedFileName


#B matrix section
BMatrixDirectory = lookupDirectory / BMatrixFileName
U1MatrixDirectory = lookupDirectory / U1MatrixFileName
U2MatrixDirectory = lookupDirectory / U2MatrixFileName
SigmaMatrixDirectory = lookupDirectory / SigmaMatrixFileName
VtMatrixDirectory = lookupDirectory / VtMatrixFileName
PseudoinverseMatrixDirectory = lookupDirectory / PseudoinverseMatrixFileName

#S section
SMatrixDirectory = lookupDirectory / SMatrixFileName

#W section
WMatrixDirectory = lookupDirectory / WMatrixFileName

#Y and Z matrix directories
YMatrixDirectory = lookupDirectory / YMatrixFileName
ZMatrixDirectory = lookupDirectory / ZMatrixFileName


#creates the class that generates the look up tables for the 
class lookUpTablesGenerator:

    #creates the initialization function
    def __init__(self,
                 M_maximum: int = numIntervalsOfInterest, #the maximum number of intervals of interest
                 highestD: int = highestDegree): #the highest degree of which to worry about

        #saves the two input variables for maximum number of intervals of interest and highest degree
        self.M_maximum = M_maximum
        self.highestD = highestD

        #we call all of the subfunction in the initialization function just to save time
        self.generateSampledBSplineLookupTables()
        self.generateIntegratedBSplineLookupTables()
        self.generateBLookupTables()
        self.generateSLookupTables()
        self.generateWLookupTables()


    #########################################################################################
    #all of these subfunctions are called in the initialization function to save the headache

    #defines the function to generate the bspline sampling functions
    def generateSampledBSplineLookupTables(self,
                                           BasisFunctionSampledFileLocation: str = BasisFunctionSampledDirectory):
        
        #instantiates the basis function writer
        sampledBasisFunctionWriter = readWriteBasisFunctions(numSamplesPerSection=numSamplesPerSection,
                                                      highestDegree=highestDegree)
        
        sampledBasisFunctionWriter.writeToNpz(BasisFunctionSampledFileLocation)

    #defines the function to generate the bspline integrated functions
    def generateIntegratedBSplineLookupTables(self,
                                              BasisFunctionIntegratedFileLocation: str = BasisFunctionIntegratedDirectory):
        #instantiates the integrator
        integratedBasisFunctionWriter = integrateBasisFunctionsContinuous(outputFileName=BasisFunctionIntegratedFileLocation,
                                                                          highestDegree = highestDegree)
        
        pass


    #creates the function to generate the B tables, and all such variations
    def generateBLookupTables(self,
                           BMatrixFileLocation: str = BMatrixDirectory,
                           U1MatrixFileLocation: str = U1MatrixDirectory,
                           U2MatrixFileLocation: str = U2MatrixDirectory,
                           SigmaMatrixFileLocation: str = SigmaMatrixDirectory,
                           VtMatrixFileLocation: str = VtMatrixDirectory,
                           PseudoinverseMatrixFileLocation: str = PseudoinverseMatrixDirectory):   


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
        np.savez(BMatrixFileLocation, **B_dictionary)
        np.savez(U1MatrixFileLocation, **U1_dictionary)
        np.savez(U2MatrixFileLocation, **U2_dictionary)
        np.savez(SigmaMatrixFileLocation, **Sigma_dictionary)
        np.savez(VtMatrixFileLocation, **Vt_dictionary)
        np.savez(PseudoinverseMatrixFileLocation, **pseudoinverseDictionary)
        
        pass
    

    #defines the function to generate the S lookup tables
    def generateSLookupTables(self,
                              SfileLocation: str = SMatrixDirectory):

        
        #creates the S dictionary
        S_dictionary = {}

        #iterates over each degree
        for d in range(self.highestD):

            #instantiates the W matrix lookUpTablesGenerator
            #the function to generate the S matrix is 
            temp_W_gen = create_W_Matrix(d=d,
                                         integratorFileName=BasisFunctionIntegratedDirectory)

            #sets the minimum M, which is always greater than the degree by at least one
            M_min = d + 1
            
            #iterates over the M's
            for M in range(M_min, (self.M_maximum + 1)):

                #gets the S matrix for the current degree and M
                current_S = temp_W_gen.S_k_M(k=d,
                                             M=M)
                
                #creates the key identifier for each current C matrix
                key = f"M{M}_d{d}"

                #stores  in the S dictionary
                S_dictionary[key] = current_S
        

        #uses savez for to save the S dictionary
        np.savez(SfileLocation, **S_dictionary)

        #finished by passing
        pass


    #creates function to generate all of the W matrices
    def generateWLookupTables(self,
                              WfileLocation: str = WMatrixDirectory):
        

        #creates the dictionary for the W data
        W_dictionary = {}

        #iterates over the valid degrees 
        for d in range(1, self.highestD + 1):

            #defines the L maximum
            L_max = d

            #iterates over the applicable L's 
            for l in range(1, (L_max + 1)):

                #defines the Minimum M
                M_min = d + 1

                #iterates over the M's
                for M in range(M_min, self.M_maximum):

                    ##creates the W temp
                    temp_W_gen = create_W_Matrix(d=d,
                                                   integratorFileName=BasisFunctionIntegratedDirectory)

                    #gets the W temp
                    W_temp = temp_W_gen.W_d_l_M(d=d,
                                                l=l,
                                                M=M)
                    
                    #creates the key
                    key = f"M{M}_d{d}_l{l}"

                    #saves to the dictionary
                    W_dictionary[key] = W_temp

        #saves the W data
        np.savez(WfileLocation, **W_dictionary)

        #returns the W dictionary
        return W_dictionary



#creates the class that reads everything
class lookUpTableReader:

    #defines the init function
    def __init__(self):

        #calls the three main functions below
        self.loadBLookupTables()
        self.loadSLookupTables()
        self.loadWLookupTables()
        pass
    


    ###########################################################################################
    #B Section of Junk
    #defines the functions to load all of the B matrices
    def loadBLookupTables(self,
                          BMatrixFile: str = BMatrixDirectory,
                          U1MatrixFile: str = U1MatrixDirectory,
                          U2MatrixFile: str = U2MatrixDirectory,
                          SigmaMatrixFile: str = SigmaMatrixDirectory,
                          VtMatrixFile: str = VtMatrixDirectory,
                          PseudoinverseMatrixFile: str = PseudoinverseMatrixDirectory):
        

    

        #loads each of those related submatrices
        B_loaded = np.load(BMatrixFile)
        U1_loaded = np.load(U1MatrixFile)
        U2_loaded = np.load(U2MatrixFile)
        Sigma_loaded = np.load(SigmaMatrixFile)
        Vt_loaded = np.load(VtMatrixFile)
        Pseudoinverse_loaded = np.load(PseudoinverseMatrixFile)

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

    #End of B Section
    ############################################################################################
    
    
    ############################################################################################
    #beginning of S Section
    #defines the function to read the S lookup npz file and load the data
    def loadSLookupTables(self,
                          SfileLocation: str = SMatrixDirectory):
        

        #loads the S matrix
        S_loaded = np.load(SfileLocation)

        #gets the S data
        self.S_data = {S_key: S_loaded[S_key] for S_key in S_loaded.files}  

        #returns that
        return self.S_data
    
    #creates the function to get the individual S
    def getIndividualS(self,
                       d: int, 
                       M: int):
        
        #creates the key for the particular S matrix
        key = f'M{M}_d{d}'


        #returns the data corresponding to the key
        return (self.S_data)[key]
    

    #end of S section
    ############################################################################################
    
    ############################################################################################
    #begin W Section

    #defines the section to load the W lookup table
    def loadWLookupTables(self,
                          WfileLocation: str = WMatrixDirectory):

        #loads the W matrix
        W_loaded = np.load(WfileLocation)

        #gets the W data
        self.W_data = {W_key: W_loaded[W_key] for W_key in W_loaded.files}

        return self.W_data
    
    #get individual W
    def getIndividualW(self,
                       M: int,
                       d: int,
                       l: int):
        
        key = f'M{M}_d{d}_l{l}'

        #returns that data
        return (self.W_data)[key]

    #end W Section
    ############################################################################################
    



#creates the class for the yz generator and reader
class YZGeneratorReader:



    #initializes
    def __init__(self):


        pass

    #creates the file to generate the tables
    def generateYZLookupTables(self,
                               YMatrixFileLocation: str = YMatrixDirectory,
                               ZMatrixFileLocation: str = ZMatrixDirectory):



        #creates the reader
        lookupReader = lookUpTableReader()

        #gets the W dictionary
        W_dict = lookupReader.W_data
        B_dict = lookupReader.B_data
        U2_dict = lookupReader.U2_data




        #creates the dictionary to store the Y data
        Y_dictionary = {}

        #creates the Z dictionary
        Z_dictionary = {}



        #iterates over all of the W matrices and corresponding 
        for W_key in W_dict:

            #calls the conversion function to ge the B key
            B_key = self.WKeyToBKey(W_key=W_key)

            #gets the corresponding U2 blocks for the B_key
            U2_temp = U2_dict[B_key]

            #gets the W_temp
            W_temp = W_dict[W_key]

            #creates the current Y matrix
            Y_temp = W_temp @ U2_temp

            #creates the current Z matrix
            Z_temp = np.transpose(U2_temp) @ W_temp @ U2_temp

            #stores the Y and Z temps in the dictionary
            #the Y and Z keys are the same as the W key, so we just reuse that
            Y_dictionary[W_key] = Y_temp
            Z_dictionary[W_key] = Z_temp

        
        #saves the Y and Z dictionaries respectively.
        np.savez(YMatrixFileLocation, **Y_dictionary)
        np.savez(ZMatrixFileLocation, **Z_dictionary)

    #defines the function to read the Y and Z dictionaries
    def loadYZTables(self,
                     YMatrixFileLocation: str = YMatrixDirectory,
                     ZMatrixFileLocation: str = ZMatrixDirectory):
        
        #loads the Y and Z matrices
        Y_loaded = np.load(YMatrixFileLocation)
        Z_loaded = np.load(ZMatrixFileLocation)


        #converts to the dictionaries
        self.Y_data = {Y_key: Y_loaded[Y_key] for Y_key in Y_loaded.files}  
        self.Z_data = {Z_key: Z_loaded[Z_key] for Z_key in Z_loaded.files}  

        return self.Y_data, self.Z_data


    #creates a function to convert a W key to a B key
    def WKeyToBKey(self, 
                W_key: str):
        
        #gets the parts of the key
        parts = W_key.split('_')

        #gets the M, d, and l values
        M = int(parts[0][1:])
        d = int(parts[1][1:])
        l = int(parts[2][1:])

        #puts them together to get the B key
        B_key = f'M{M}_d{d}'

        return B_key