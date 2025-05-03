#this file implements the efficient implementation of the other functions in the legacy matrix helpers file

import numpy as np

#imports the function to evaluate the uniform basis function at a specific set of points
from eVTOL_BSplines.path_generation_helpers.matrix_helpers import uniform_basis_function_evaluation, D_d_M, D_d_l_M
import os, sys
from pathlib import Path
from scipy.integrate import quad


#creates the class to sample the basis functions of the 
class basisFunctionSampler:

    #creates the initialization
    def __init__(self,
                 degree: int, #saves the degree of this thing
                 numSegmentsPerSection: int): #saves the number of samples per section of the thing
        
        #seaves the degree for this thing
        self.degree = degree

        #saves the number of sections, which is one plus the degree
        self.numSections = self.degree + 1

        #saves the number of samples per section, which is one plus the number of segments. (we want to get the ends.)
        self.numSamplesPerSection = numSegmentsPerSection

        #gets the Ts
        self.Ts = (1.0/float(numSegmentsPerSection))
        

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

                if currentTime > 0.99:
                    

                    vegetable = 0

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
                                                  numSegmentsPerSection=self.numSamplesPerSection)
            
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

#creates the class to read from an npz file
class readBasisFunctions:

    def __init__(self, 
                 fileName: str,
                 highestDegree: int):
        
        #loads the npz
        loadedFile = np.load(fileName)

        #saves the highest degree
        self.highestDegree = highestDegree

        self.numBasisTypes = self.highestDegree + 1

        #saves the 
        self.loaded_list = [loadedFile[f"array_{i}"] for i in range(self.numBasisTypes)]

    #returns the loaded list
    def getLoadedList(self):
        return self.loaded_list

#'''
#creates the class to integrate the subsections of the basis functions
class integrateBasisFunctionsSampled:

    def __init__(self,
                 fileName: str,
                 highestDegree: int):

        self.highestDegree = highestDegree


        #loads in the file which has all of the sections 
        self.fileName = fileName

        temp1 = os.fspath(Path(__file__).parents[0])
        temp2 = os.path.abspath(os.path.join(temp1, fileName))

        #creates the reader/writer
        reader = readBasisFunctions(fileName=self.fileName, highestDegree=highestDegree)

        #gets the loadedList
        self.loadedList = reader.getLoadedList()

        
    #creates the function to iterate through the list and perform all of the integreations

    def createIntegrations(self,
                           outputFileName: str):#saves the outputFile name


        #creates the list to store the integrated parts for each degree
        degreeIntegratedParts = []

        #iterates over the degrees
        for i in range(self.highestDegree + 1):
            

            #gets the item in the list
            sectionList = (self.loadedList)[i]

            #gets the number of subsections (which is degree plus 1)
            numSubsections = i + 1

            #sets the number of integrations as the number of subsections squared
            numIntegrations = numSubsections**2
            #the number of combinations is the number of sections squared
            

            #gets the 2d integrated list
            wholeIntegratedList = []


            #next we need to work on adding integration
            #iterates over the number of subsections, where this number refers to the first one (thet static one)
            for j in range(numSubsections):
                
                currentIntegratedList = []
                #iterates over the moving subsections
                for k in range(numSubsections):

                    #gets the first subsections
                    subsection_1 = sectionList[i]
                    #gets the second subsection
                    subsection_2 = sectionList[j]

                    #gets the Ts value
                    Ts = 1.0/(np.size(subsection_1))

                    #gets the multiplication of subsection 1 and 2
                    subsectionProduct = subsection_1*subsection_2

                    #gets the trapezoidal integration product
                    integral = np.trapz(y=subsectionProduct, dx=Ts)
                    
                    #appends the integral
                    currentIntegratedList.append(integral)
                #appends the currentIntegratedList
                wholeIntegratedList.append(currentIntegratedList)

            #gets the array equivalent of the whole integrated list
            wholeIntegratedArray = np.array(wholeIntegratedList)
        
            #appends the whole integrated list (for the current degree polynomial) to the integrated parts list
            degreeIntegratedParts.append(wholeIntegratedArray)

            #iterates over each number of subsections


        #creates the dictionary of the integration arrays
        arrays_dictionary = {f'degree_{i}': array for i, array in enumerate(degreeIntegratedParts)}

        #saves to an npz file
        np.savez(outputFileName, **arrays_dictionary)

        quesadilla = 0
#'''

class integrateBasisFunctionsContinuous:

    def __init__(self,
                 outputFileName: str,
                 highestDegree: int):
        

        self.highestDegree = highestDegree
        self.outputFileName = outputFileName

        temp1 = os.fspath(Path(__file__).parents[0])
        temp2 = os.path.abspath(os.path.join(temp1, outputFileName))


        #calls the create integrations function
        self.integrations_dict = self.createIntegrations()

        np.savez(outputFileName, **(self.integrations_dict))

    #creates the function to get the integrations
    def createIntegrations(self):

        #creates the function to get the degree Integrated Parts
        integratedSectionsByDegree = []

        #iterates over each degree 
        for currentDegree in range(self.highestDegree + 1):

            
            #creates the individual degree integrated list
            currentDegreeList = []
            

            #the number of subsections of the function is i + 1
            numSubsections = int(currentDegree+1)

            #sets the number of subsections to iterate over
            for j in range(numSubsections):


                #gets the currenttemp list
                currentTempList = []
                #iterates over the moving subsection
                for k in range(numSubsections):
                    #gets the time offset for the first function
                    offset_1 = float(j)
                    offset_2 = float(k)
                    #gets the current result and erro
                    integralResult, error = quad(self.multipliedFunction, 
                                                 a=0,
                                                 b=1,
                                                 args=(currentDegree,offset_1,offset_2))
                    
                    #appends the integral REsult to the list
                    currentTempList.append(integralResult)
                currentDegreeList.append(currentTempList)
            

            currentDegreeArray = np.array(currentDegreeList)
            #appends the current Degree List to the full integrated List
            integratedSectionsByDegree.append(currentDegreeArray)

        #creates the dictionary
        arrays_dictionary = {f'degree_{i}': array for i, array in enumerate(integratedSectionsByDegree)}

        return arrays_dictionary


    


    #creates the function definition for the same function and two different times multiplied together
    def multipliedFunction(self, 
                           time: float,
                           degree: int,
                           offset_1: float,
                           offset_2: float):
        
        time_1 = time + offset_1
        time_2 = time + offset_2
        
        value_1 = uniform_basis_function_evaluation(time=time_1, degree=degree, alpha=1.0)
        value_2 = uniform_basis_function_evaluation(time=time_2, degree=degree, alpha=1.0)

        #gets the multiplication product between the two different ones
        product = value_1*value_2

        #returns the product
        return product







#creates the class to efficiently construct the W matrix
#this class needs to Create the S_k_M matrix. 
class create_W_Matrix:

    #creates the init function
    def __init__(self,
                 d: int, #degree of the basis splines
                 integratorFileName: str):
        
        temp1 = os.fspath(Path(__file__).parents[0])
        temp2 = os.path.abspath(os.path.join(temp1, integratorFileName))

        self.numBasisTypes = 6
         

        #saves the d degree
        self.d = d

        #reads in the file
        loadedFile = np.load(temp2) 
        #saves the 
        self.integration_Matrix_List = [loadedFile[f"degree_{i}"] for i in range(self.numBasisTypes)]

        #gets the sum along the diagonals for the first column
        columnDiagonalSum = []


        
        #gets the sum along the diagonals for the first row
        rowDiagonalSum = []




        chepeta = 0


    def S_k_M(self,
              k: int,
              M: int):
        #creates the function to create the S_k_M matrix
        

        #gets the correct integration matrix
        integrationMatrix = (self.integration_Matrix_List)[k]


        #creates a test matrix to test whether we are getting the correct number of sections
        numSectionsTestMatrix = np.zeros((k+M, k+M))

        length = k + M
        #creates the matrix
        S_matrix = np.zeros((length, length))
        
        #iterates through for the whole thing.
        for i in range(length):
            for j in range(length):


                #gets the difference between i and j
                difference = i-j


                #checks if the difference is greater than the degree
                if np.abs(difference) > k:
                    S_matrix[i,j] = 0.0

                    numSectionsTestMatrix[i,j] = 0

                
                #case it is in the cutoff area
                elif (i < k and j < k) or (i >= M and j >= M):
                    
                    #gets the number of sections of removal from this from 
                    numRemovals = self.getNumRemovals(M=M,
                                                      k=k,
                                                      i=i,
                                                      j=j)
                    

                    #sets the difference index
                    differenceIndex = np.abs(difference)

                    #gets the sign of the difference index
                    if np.sign(difference) == 1:
                        isColumn = True
                    else:
                        isColumn = False

                    #sets the index for the diagonal section as 
                    #gets the full diagonal section
                    fullDiagonalSection = self.getIMatrixDiagonal(I=integrationMatrix,
                                                                  index=differenceIndex,
                                                                  isColumn=isColumn) 


                    if (i < k and j < k):
                        removingLeft = True
                    elif (i >= M and j >= M):
                        removingLeft = False
                    
                    #gets the section of the diagonal for this one, 
                    # and this one is taking off the left side
                    partialSection = self.getDiagonalSection(removingLeft=removingLeft,
                                                             numRemovals=numRemovals,
                                                             diagonal=fullDiagonalSection)
                    
                    #gets the sum of the partial Section Sum
                    partialSectionSum = np.sum(partialSection)

                    partialSectionLength = np.size(partialSection)

                    #sets the portions of the matrix
                    S_matrix[i,j] = partialSectionSum

                    numSectionsTestMatrix[i,j] = partialSectionLength
                
                #otherwise, we get the whole diagonal and sum it up
                else:
                    #sets the difference index
                    differenceIndex = np.abs(difference)

                    #gets the sign of the difference index
                    if np.sign(difference) == 1:
                        isColumn = True
                    else:
                        isColumn = False

                    #sets the index for the diagonal section as 
                    #gets the full diagonal section
                    fullDiagonalSection = self.getIMatrixDiagonal(I=integrationMatrix,
                                                                  index=differenceIndex,
                                                                  isColumn=isColumn) 
                    

                    #gets the full sum and the full sum length
                    fullSum = np.sum(fullDiagonalSection)
                    fullLength = np.size(fullDiagonalSection)

                    #saves them
                    S_matrix[i,j] = fullSum
                    numSectionsTestMatrix[i,j] = fullLength


        #returns the S matrix
        return S_matrix                   



    #defines the function to sum along a particular diagonal corresponding to starting on the top row
    def sumDiagonalTopRow(self,
                          A: np.ndarray, #the A matrix on which to perform the diagonal summation
                          index: int):#the starting index to work on


        #gets the diagonal length
        diagonalLength = self.getDiagonalLength(A=A, index=index)

        #creates the sum as zero
        sum = 0.0

        #iterates over all the elements along the diagonal
        for i in range(diagonalLength):

            #gets the current value
            currentValue = A[i,i+index]
            #adds to the sum
            sum = sum + currentValue


        #returns the summation
        return sum
    
    #does the same sum, but for the left column
    def sumDiagonalLeftConlum(self,
                              A: np.ndarray,
                              index: int):
        #gets the diagonal length
        diagonalLength = self.getDiagonalLength(A=A, index=index)
        
        #initializes the sum
        sum = 0.0
        #iterates through and sums
        for i in range(diagonalLength):

            currentValue = A[i+index, i]
            sum = sum + currentValue

        #returns the summation
        return sum


    #defines the function to get the length of a particular diagonal
    def getDiagonalLength(self,
                          A: np.ndarray,
                          index: int): #index is the offset off the diagonal

        #gets the shape of the matrix A
        A_shape = np.shape(A)

        #checks if we have the incorrect shape of a matrix
        if A_shape[0] != A_shape[1]:

            raise ValueError("Matrix is not square")
        
        #checks if we have a valid index
        if index >= A_shape[0]:

            raise ValueError("Index is larger than matrix")
        
        if index < 0:
            raise ValueError("Index Cannot be Less than One")


        #the return value is the length of A_shape[0] minus index
        diagonalLength = A_shape[0] - index

        return diagonalLength
    
    #defines the function to get the diagonal of the I matrix as an array (flattened, of course)

    def getIMatrixDiagonal(self,
                           I: np.ndarray, #the I matrix in question
                           index: int, #the starting index
                           isColumn: bool = True):#bool that defines whether the index corresponds to a column (if false, then we're getting for a row) this is important if it's not symmetric

        #gets the diagonal length
        diagonalLength = self.getDiagonalLength(A=I, index=index)
        
        section = []

        #gets the elements of the diagonal section
        for i in range(diagonalLength):
            
            #case is Column
            if isColumn:
                currentValue = I[i+index,i]
            #case is not column
            else:
                currentValue = I[i, i+index]

            #adds the current Vlue to the 
            section.append(currentValue)

        #converts to the section
        section = np.array(section)

        return section
        

    #creates the function to get a specified section of the diagon
    def getDiagonalSection(self, 
                           removingLeft: bool, #boolean to indicate whether we are removing from the left side, or alternativelym from the right side
                           numRemovals: int, #the number of removals from the 
                           diagonal: np.ndarray):
        

        #gets the length of the diagonal
        length = np.size(diagonal)

        #if the number of removals is greater than or equal to the length raise an error
        #or if it's negative
        if numRemovals >= length or numRemovals < 0:
            raise ValueError(f"Num Removals must be greater that 0, and less than {length}")

        #gets the section
        #case removing from the left side

        if removingLeft:

            #we take off from the left side the number to be removed
            output = diagonal[numRemovals:]

        #otherwise, takeoff from the end
        else:
            output = diagonal[:(length-numRemovals)]

        #returns the output
        return output
    

    #creates a function to determine the number of removals 
    def getNumRemovals(self,
                       M: int, #the number of intervals 
                       k: int, #the degree variable
                       i: int, #iteration variable
                       j: int): #iteration variable
        
        #if j is less than or equal to i

        #case i is less than or equal to d
        if i <= k and j <= k:
            #this block is the necessary conditions to determine what the number of sections to remove is
            if j <= i:
                numRemovals = k - i
            else:
                numRemovals = k - j

        elif i >= M and j >= M:
            #sets the offset so that the ending block looks like the beginning block above. It makes the logic more logical
            offset = M + k - 1
            #gets the new i and j
            i_new = offset - i
            j_new = offset - j

            #this block is the necessary conditions to determine what the number of sections to remove is
            if j_new <= i_new:
                numRemovals = k - i_new
            else:
                numRemovals = k - j_new

        #returns the number of removals
        return numRemovals



        
