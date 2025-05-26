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
                           pseudoinverseLocation: str = "lookUpTables/PseudoinverseMatrices.npz",
                           fileLocation: str = "lookUpTables/B_Matrices.npz"):

        temp1 = os.fspath(Path(__file__).parents[0])
        intputFilePath = os.path.abspath(os.path.join(temp1, fileLocation))

        pseudoinverseFilePath = os.path.abspath(os.path.join(temp1, pseudoinverseLocation))

        #creates the B lookup tables, with the svd data
        B_dict = {}

        #creates the B metadata
        B_metadata = {}

        #creates the pseudoinverse lookup table
        pseudoinverse_dict = {}
        
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


                #appends the key to the metadata
                B_metadata[mainKey] = np.array([degree, M])

                #obtains the pseudoinverse and saves it with the key in the pseudoinverse dictionary
                V = np.transpose(Vt)
                
                Sigma_inv_items = []
                #gets the reciprocal of Each item in S
                for item in S:
                    temp = 1.0/item
                    #saves this to the items list
                    Sigma_inv_items.append(temp)
                
                #turns the Sigma inverse items into a diagonal matrix
                Sigma_inv = np.diag(np.array(Sigma_inv_items))

                #gets the U1 transpose
                U1_T = np.transpose(U1)

                #gets the pseudoinverse of B
                B_cat_pseudoinverse = V @ Sigma_inv @ U1_T
                #appends it to the pseudoinverse list
                pseudoinverse_dict[mainKey] = (B_cat_pseudoinverse)



        #saves the file to that file path
        np.savez(intputFilePath, **B_dict, metadata=B_metadata)
        np.savez(pseudoinverseFilePath, **pseudoinverse_dict)

        return B_dict, B_metadata


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
                              fileLocation: str = "lookUpTables/W_d_l_M_Matrices.npz",
                              W_metadataLocation: str = "lookUpTables/W_metadata.npz"):

        temp1 = os.fspath(Path(__file__).parents[0])
        inputFilePath = os.path.abspath(os.path.join(temp1, fileLocation))
        
        metadataFilePath = os.path.abspath(os.path.join(temp1, W_metadataLocation))


        #crates the dictionary to store all of the W matrices
        W_dict = {}

        W_metadata = {}
        
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

                    W_metadata[key] = np.array([degree, ell, M])



        #saves the W_dictionary unpacked to an npz
        np.savez(inputFilePath, **W_dict)
        np.savez(metadataFilePath, **W_metadata)
        
        #passes and does not return anything
        return W_dict, W_metadata

#''''''
    def generateYZLookupTables(self,
                               y_fileLocation: str = "lookUpTables/Y_Matrices.npz",
                               z_fileLocation: str = "lookUpTables/Z_Matrices.npz"):

        temp1 = os.fspath(Path(__file__).parents[0])
        y_inputFilePath = os.path.abspath(os.path.join(temp1, y_fileLocation))
        z_inputFilePath = os.path.abspath(os.path.join(temp1, z_fileLocation))    


        #creates the reader
        lookupReader = lookUpTableReader()
        #gets the W dictionary
        W_dict = lookupReader.W_data
        #gets the B dictionary
        B_dict = lookupReader.B_data
        #gets the metadata key for everything
        metadata = lookupReader.W_metadata

        #creates the dictionary for the Y's
        #remember that Y_M_d_l = W_M_d_l * U2_M_d
        Y_dict = {}

        #creates the dictionary for the Z's
        #remember that Z_M_d_l = U2_M_d^T * W_M_d_l * U2_M_d
        Z_dict = {}

        #iterates over each element in the W dictionary
        for d_l_M_value in metadata.values():
            
            d = d_l_M_value[0]
            ell = d_l_M_value[1]
            M = d_l_M_value[2]
            #recreates the key for the W
            W_key = f'degree{d}_l{ell}_M{M}'
            
            #recreates the B main key
            B_main_key = f"d{d}_M{M}"

            #recreates the key for U2
            U2_key = f"{B_main_key}_U2"

            #gets the W matrix from the W Key
            W_temp = W_dict[W_key]
            #gets the U2 Matrix from the 
            U2_temp = B_dict[U2_key]


            #creates the Y_temp
            Y_temp = W_temp @ U2_temp

            #creates the Z_temp to store in the dictionary
            Z_temp = U2_temp.T @ W_temp @ U2_temp

            #stores them in the dictionaries
            Y_dict[W_key] = Y_temp
            Z_dict[W_key] = Z_temp


            potato = 0

        #saves the Y_dict and Z_dicts to the respective arrays
        np.savez(y_inputFilePath, **Y_dict)
        np.savez(z_inputFilePath, **Z_dict)


        return Y_dict, Z_dict
    
    def generatePseudoinverseLookupTables(self,
                                          fileLocation: str = "lookUpTables/Pseudoinverses.npz"):
        potato = 0

#'''




#creates the class to read from and process the look up table files
class lookUpTableReader:

    #creates the init function
    def __init__(self):
        #calls the three load functions to start by loading them up
        self.loadWLookupTable()
        self.loadSLookupTables()
        self.loadBLookupTable()
        self.loadWMetadata()
        self.loadYLookupTable()
        self.loadZLookupTable()
        self.loadPseudoinverseLookupTable()

        pass

    #creates the function to read the S Lookup tables
    def loadSLookupTables(self,
                          fileLocation: str = "lookUpTables/S_Matrices.npz"):
    
        temp1 = os.fspath(Path(__file__).parents[0])
        inputFilePath = os.path.abspath(os.path.join(temp1, fileLocation))

        S_loaded = np.load(inputFilePath)

        #gets the loaded S matrix
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
    def loadBLookupTable(self,
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
    
    #creates the function to load the pseudoinverse lookup table
    def loadPseudoinverseLookupTable(self,
                                     fileLocation: str = 'lookUpTables/PseudoinverseMatrices.npz'):
        
        temp1 = os.fspath(Path(__file__).parents[0])
        inputFilePath = os.path.abspath(os.path.join(temp1, fileLocation))

        #gets the loaded pseudoinverse array
        PseudoinverseLoaded = np.load(inputFilePath)

        #gets the data
        PseudoinverseData = {key: PseudoinverseLoaded[key] for key in PseudoinverseLoaded.files}

        #saves the data
        self.PseudoinverseData = PseudoinverseData

        #returns the data
        return PseudoinverseData
    
    #creates the function to get an individual bit of pseudoinverse data
    def getIndividualPseudoinverse(self, d: int, M: int):

        #creates the key
        key = f'd{d}_M{M}'

        #gets the B_inv
        B_inv = (self.PseudoinverseData)[key]

        #returns that
        return B_inv

    
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
    def loadWLookupTable(self,
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
    
    #function to moad the metadata file
    def loadWMetadata(self,
                      fileLocation: str = "lookUpTables/W_metadata.npz"):
        
        temp1 = os.fspath(Path(__file__).parents[0])
        inputFilePath = os.path.abspath(os.path.join(temp1, fileLocation))
        
        W_meta_loaded = np.load(inputFilePath)

        W_metadata = {key: W_meta_loaded[key] for key in W_meta_loaded.files}

        self.W_metadata = W_metadata

        return W_metadata
    
    def getIndividualWMetadata(self, key: str):

        #gets the associated array for the input key
        array = (self.W_metadata)[key]

        degree = array[0]
        ell = array[1]
        M = array[2]

        return degree, ell, M


    #defines the function to load the Y matrices
    def loadYLookupTable(self,
                         fileLocation: str = "lookUpTables/Y_Matrices.npz"):
        temp1 = os.fspath(Path(__file__).parents[0])
        y_inputFilePath = os.path.abspath(os.path.join(temp1, fileLocation))

        #loads it from that file path
        Y_data_loaded = np.load(y_inputFilePath)

        #converts it back to a dictionary
        Y_data = {key: Y_data_loaded[key] for key in Y_data_loaded.files}

        #saves it to self.
        self.Y_data = Y_data

        return Y_data
    
    #defines the same for the Z matrices
    def loadZLookupTable(self,
                         fileLocation: str = "lookUpTables/Z_Matrices.npz"):
        temp1 = os.fspath(Path(__file__).parents[0])
        z_inputFilePath = os.path.abspath(os.path.join(temp1, fileLocation))

        #loads it from that file path
        Z_data_loaded = np.load(z_inputFilePath)

        #converts it back to a dictionary
        Z_data = {key: Z_data_loaded[key] for key in Z_data_loaded.files}

        #saves it to self.
        self.Z_data = Z_data

        return Z_data