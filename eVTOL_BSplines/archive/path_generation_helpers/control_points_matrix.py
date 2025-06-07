#creates the third attempt at an efficient
#matrix generator for the control points


import numpy as np
import os, sys
from pathlib import Path
#imports the things I have 
from eVTOL_BSplines.path_generation_helpers.matrix_helpers import D_d_l_M, B_init_final, B_init_final_svd




class controlPointsGenerator:

    #creates the initialization function
    def __init__(self,
                 dimension: int,
                 degree: int,
                 M: int
                 ):

        #initializes the create W matrix class
        self.W_matrix_creator = create_W_Matrix(d=degree,
                                                integratorFileName= "lookUpTables/degree_5_integrations.npz")
        
        #creates the B init and ifnal
        B = B_init_final(degree=degree, alpha=1.0, M=M)

        self.degree = degree
        self.M = M



    #creates the function to generate the control points
    def generateControlPoints(self,
                              Start: np.ndarray,
                              End: np.ndarray,
                              rho: np.ndarray):
        
        #creates the matrix of the boundary conditions
        boundaryConditions = np.concatenate((Start, End), axis=1)
        
        #gets the W matrix
        W_d_M_rho = (self.W_matrix_creator).W_d_M_rho(d=self.degree,
                                                      M=self.M,
                                                      rho=rho)
        
        #gets the B matrix
        U1, U2, Sigma, V_T = B_init_final_svd(degree=self.degree,
                                              alpha=1.0,
                                              M=self.M)
        
        U1_T = np.transpose(U1)
        U2_T = np.transpose(U2)
        Sigma_inv = np.linalg.inv(Sigma)
        V = np.transpose(V_T)
        
        #this section solves for the Z matrix
        #gets the A1 matrix
        A1 = np.transpose(U2) @ W_d_M_rho @ U2

        #gets the B1 matrix
        B1 = - boundaryConditions @ V @ Sigma_inv @ U1_T @ W_d_M_rho @ U2

        #gets the transposed verstion of each of those matrices
        A1_T = np.transpose(A1)
        B1_T = np.transpose(B1)

        #uses the lstsq function to get the Z_T solution
        Z_T, residuals, rank, s = np.linalg.lstsq(A1_T, B1_T, rcond=None)

        #gets the Z matrix
        Z = np.transpose(Z_T)


        #obtains the Control Points
        C = boundaryConditions @ V @ Sigma_inv @ U1_T + Z @ U2_T


        #returns the Control points
        return C



    
#creates the class to efficiently construct the W matrix
#this class needs to Create the S_k_M matrix. 
class create_W_Matrix:

    #creates the init function
    def __init__(self,
                 d: int): #degree of the basis splines
                 
        
        temp1 = os.fspath(Path(__file__).parents[0])
        intputFilePath = os.path.abspath(os.path.join(temp1, 'lookUpTables/degree_5_integrations.npz'))

        self.numBasisTypes = 6
         

        #saves the d degree
        self.d = d

        #reads in the file
        loadedFile = np.load(intputFilePath) 
        #saves the 
        self.integration_Matrix_List = [loadedFile[f"degree_{i}"] for i in range(self.numBasisTypes)]

        chepeta = 0

    #defines the function to obtain the matrix S
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

                    #saves them
                    S_matrix[i,j] = fullSum


        #returns the S matrix
        return S_matrix


    #defines the function to obtain the W_d_l_M matrix (which is the individual matrix
    def W_d_l_M(self,
                d: int, #defines the degree of the bspline
                l: int, #defines the degree of the derivative
                M: int): #defines the number of intervals of interest
        

        #gets the D matrix
        D = D_d_l_M(d=d,l=(l-1),M=M)
        
        #gets the S matrix
        S = self.S_k_M(k=(d-l),M=M)

        W = D @ S @ np.transpose(D)

        #returns the W matrix
        return W


    #defines the function to get the proper subsections of the matrices
    def W_applicable_subsections(self,
                                 d: int,
                                 M: int,
                                 rho: np.ndarray):
        
        W_partition = self.W_partitioned(d=d, M=M, rho=rho)

        #gets the applicable matrices
        W_cent_cent = (W_partition[1])[1]

        #gets the matrix top center
        W_top_cent = (W_partition[0])[1]

        #gets the matrix right center
        W_cent_right = (W_partition[1])[2]

        #returns these
        return W_cent_cent, W_top_cent, W_cent_right

    #defines the function to obtain the W_d_M_rho function,
    #but it obtains it as a partitioned matrix
    def W_partitioned(self,
                      d: int, #the degree of the polynomial
                      M: int, #the number of intervals of interest
                      rho: np.ndarray): #the scaling summation array

        #gets the W unpartitioned matrix from the other function
        W_unpartitioned = self.W_d_M_rho(d=d,
                                         M=M,
                                         rho=rho)


        #gets the partitioned matrix of W
        W = self.getWMatrixPartition(W=W_unpartitioned,
                                     d=d,
                                     M=M)

        #returns the W matrix
        return W

                

    #defines the function to obtain the complete W_d_M_rho matrix
    def W_d_M_rho(self,
                  d: int, #the degree of the bspline to be created
                  M: int, #the number of intervals of interest
                  rho: np.ndarray): #the scaling matrix used to scale the various derivatives
        
        #iterates by the rho length
        L = np.size(rho)

        #getsthe sum matrix
        sumMatrix = np.zeros((d+M, d+M))

        for l in range(L):
            #gets the scaled W
            scaledW = rho.item(l) * self.W_d_l_M(d=d,l=l,M=M)


            sumMatrix = sumMatrix + scaledW

        #returns the Sum Matrix

        return sumMatrix

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
    
    #gets the partition of a matrix based on d and M
    def getWMatrixPartition(self,
                           W: np.ndarray, #the main array from which to get the partition
                           d: int, #d: the degree of the bspline associated with the Matrix
                           M: int): #M: the number of intervals of interest
        
        #checks if this is the correct shape of matrix:
        W_shape = np.shape(W)

        desiredShape = (d+M, d+M)

        if W_shape != desiredShape:

            raise ValueError(f"Shape of W: {W_shape} is incompatible with specified shape: {desiredShape}")

        #obtains the nine subsections of the matrix
        #left column of sections
        W_00 = W[:d,:d]
        W_10 = W[d:M,:d]
        W_20 = W[M:,:d]
        #middle column of sections
        W_01 = W[:d,d:M]
        W_11 = W[d:M,d:M]
        W_21 = W[M:,d:M]
        #right column of sections
        W_02 = W[:d,M:]
        W_12 = W[d:M,M:]
        W_22 = W[M:,M:]

        left = [W_00, W_10, W_20]
        center = [W_01, W_11, W_21]
        right = [W_02, W_12, W_22]

        W_partitioned = [left, center, right]

        return W_partitioned
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







