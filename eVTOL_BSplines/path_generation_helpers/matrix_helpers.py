#this file is for helping out with the matrix representation. 
#I wish that some of the functionality had been included in the path generator files. 
#alas, they were not, but that's okay for now.

import os, sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from bsplinegenerator.bsplines import BsplineEvaluation
import math
from typing import List

import scipy.integrate as integrate


from bsplinegenerator.table_evaluation import cox_de_boor_table_basis_function

sys.path.insert(0,os.fspath(Path(__file__).parents[2]))

temp1 = os.fspath(Path(__file__).parents[0])
temp2 = os.path.abspath(os.path.join(temp1, 'submodules/path_generator'))
sys.path.insert(0,temp2)
tempPath = sys.path

from eVTOL_BSplines.submodules.path_generator.path_generation.matrix_evaluation import *
from eVTOL_BSplines.submodules.path_generator.path_generation.matrix_evaluation import *
from eVTOL_BSplines.submodules.path_generator.path_generation.obstacle import *
from eVTOL_BSplines.submodules.path_generator.path_generation.path_generator import *
from eVTOL_BSplines.submodules.path_generator.path_generation.path_plotter import *
from eVTOL_BSplines.submodules.path_generator.path_generation.safe_flight_corridor import *
from eVTOL_BSplines.submodules.path_generator.path_generation.waypoint_data import *


from bsplinegenerator.table_evaluation import cox_de_boor_table_basis_function, __cox_de_boor_table_basis_function_whole_table



'''
    def __create_knot_points(self):
        
        This function creates evenly distributed knot points
        
        number_of_control_points = self._num_control_points
        number_of_knot_points = number_of_control_points + self._order + 1
        knot_points = np.arange(number_of_knot_points)*self._scale_factor + self._start_time - self._order*self._scale_factor
        return knot_points
#'''

#creates the uniform knot point generator
def uniform_knot_point_generator(M: int, #number of intervals of interest
                                 degree: int, #the degree of the bspline functions
                                 alpha: float, #the factor of scaling for the knot vector
                                 start_time: float): #the start time for the zeroth knot point
    
    #gets the number of knot points
    num_knot_points = M + 2*degree + 1

    #gets the knot points, which is in the range of zero to the number of knot points
    #then it is shifted back by the degree
    #and then it is scaled by the scaling factor
    #finally, the start time defines the initial start time shift
    knot_points = (np.arange(num_knot_points) - degree)*alpha + start_time


    #returns the knot_points
    return knot_points




#creates the function to get D matrices
def D_d_M(d: int,  #the degree of the spline
              M: int): #the number of intervals of interest

    #creates the first zeros
    zeros = np.zeros((1,M+d-1))

    identityMatrix = np.eye(M+d-1)
    #gets the first matrix
    A_1 = np.concatenate((zeros, identityMatrix), axis=0)

    #creates the second matrix
    A_2 = np.concatenate((identityMatrix, zeros),axis=0)

    #gets the D matrix
    D_d_M = A_1 - A_2
    #returns it
    return D_d_M

#creates the uniform cox_de_boor basis function
def uniform_basis_function_evaluation(time: float, #the time at which to evaluate the function
                                      degree: int, #the degree of the polunomial we are using
                                      alpha: float): #the scaling factor alpha, which in this case just scales time
    
    #gets the scaled time
    time_scaled = time/alpha

    #if the input degree is zero
    if degree == 0:
        #case the input time is between 0 and 1
        if 0.0 <= time_scaled and time_scaled < 1.0:
            basis_function = 1.0
        #case the input time is outside of that bound
        else:
            basis_function = 0.0
    #otherwise, we go into a recursion
    else:
        #first term section
        term_1 = (time_scaled/degree)*uniform_basis_function_evaluation(time=time_scaled,
                                                                  degree=(degree-1),
                                                                  alpha=alpha)
        #second term section
        term_2 = ((degree+1-time_scaled)/(degree))*uniform_basis_function_evaluation(time=(time_scaled-1),
                                                                               degree=(degree-1),
                                                                               alpha=alpha)
        
        #sets the basis function to the sum of the two terms
        basis_function = term_1 + term_2

    #returns the basis function
    return basis_function


#creates the function to get the b_d_M(t)

def b_d_M_t_vector(time: float, #the current evaluation time
                   degree: int, #the degree of the spline function
                   alpha: float, #the scaling function
                   M: int): #the number of intervals of interest in the function
    
    #gets the shape of the array
    length = degree + M
    #creates the matrix to store the values
    b_d_M = np.zeros((length,1))

    #iterates through the length and gets the evaluation using the uniform function evaluation
    for i in range(length):
        #creates the input time
        input_time = time + degree - i
        #gets the value from the function evaluation
        value = uniform_basis_function_evaluation(time=input_time,
                                                  degree=degree,
                                                  alpha=alpha)
        
        #puts the value in its place
        b_d_M[i,0] = value
    

    #returns the matrix
    return b_d_M



#creates the function to create the B matrix, as well as the B_hat matrix
def B_M_matrix(time: float, #time at which to evaluate the matrix
                   degree: int, #the degree at which to evaluate the matrix
                   alpha: float, #the scaling factor
                   M: int): #the number of intervals of interest
    #creates the matrix itself
    B_d_M = np.zeros((M+degree, degree))

    #iterates through each column
    for i in range(degree):
        #gets the degree for the b_d_M vector
        b_degree = (degree-i)

        #gets the applicable b vector
        b_d_M = b_d_M_t_vector(time=time,
                               degree=b_degree,
                               alpha=alpha,
                               M=M)
        
        #sets the initial net D matrix, which we create as an identity matrix
        D = np.eye(degree + M)
        #iterates through to get the net D matrix
        for j in range(i):
            #gets the D degree
            D_degree = degree - j
            #gets the next temp d matrix
            D_temp = D_d_M(d=D_degree,
                               M=M)
            
            #multiplies the D_temp into the whole D Matrix
            D = D @ D_temp
        
        #now that we have the D matrix, we multiply it to the b vector
        b_vector_result = D @ b_d_M

        b_vector_result = b_vector_result.reshape((np.size(b_vector_result)))

        #puts it in the correct slot
        B_d_M[:,i] = b_vector_result

    ##########################################################################
    #section of the code that gets the B_hat vector

    #converts the time to an int and uses that as an index
    time_index = int(time)
    #then gets the appropriate section.
    B_hat_d_M = B_d_M[(time_index):(time_index+degree),:]
        
    #returns the whole matrix, and the B_hat matrix
    return B_d_M


#creates function that creates the B_hat and B_hat inverse matrix for a particular thing
def B_hat_B_hat_inv(degree: int, #the degree at which to evaluate the matrix
                   alpha: float, #the scaling factor
                   M: int): #the number of intervals of interest
    
    #sets the time time
    time = 0.0
    
    #creates the matrix itself
    B_d_M = np.zeros((M+degree, degree))

    #iterates through each column
    for i in range(degree):
        #gets the degree for the b_d_M vector
        b_degree = (degree-i)

        #gets the applicable b vector
        b_d_M = b_d_M_t_vector(time=time,
                               degree=b_degree,
                               alpha=alpha,
                               M=M)
        
        #sets the initial net D matrix, which we create as an identity matrix
        D = np.eye(degree + M)
        #iterates through to get the net D matrix
        for j in range(i):
            #gets the D degree
            D_degree = degree - j
            #gets the next temp d matrix
            D_temp = D_d_M(d=D_degree,
                           M=M)
            
            #multiplies the D_temp into the whole D Matrix
            D = D @ D_temp
        
        #now that we have the D matrix, we multiply it to the b vector
        b_vector_result = D @ b_d_M

        b_vector_result = b_vector_result.reshape((np.size(b_vector_result)))

        #puts it in the correct slot
        B_d_M[:,i] = b_vector_result

    ##########################################################################
    #section of the code that gets the B_hat vector

    #converts the time to an int and uses that as an index
    time_index = int(time)
    #then gets the appropriate section.
    B_hat_d_M = B_d_M[(time_index):(time_index+degree),:]


    #gets the inverse matrix
    B_hat_d_M_inv = np.linalg.inv(B_hat_d_M)
        
    #returns the whole matrix, and the B_hat matrix
    return B_hat_d_M, B_hat_d_M_inv


###################################################################################
#this is the section to get the B hat matrix, but in the simplified manner, and not the complex one


#gets the B_d Vector
def b_d_vector(degree: int,
               alpha: float = 1.0):

    #in this case, length = degree
    length = degree
    #greates the matrix to store the values
    b_d = np.zeros((length, 1))

    #iterates through the length and gets the evaluatio nusing the uniform function evaluation
    for i in range(length):
        ##creates the input time from top to bottom of matrix
        input_time = degree - i

        #gets  the evaluation  value at that point
        value = uniform_basis_function_evaluation(time=input_time,
                                                  degree=degree,
                                                  alpha=alpha)
        
        #puts the value in its place
        b_d[i,0] = value

    #returns the vector
    return b_d


#creates a function to get a 

#creates the simplified function to obtain just the B_hat matrix, 
#which is constant, regardless of the time
def B_hat_B_hat_inv_simplified(degree: int,
                               alpha: float):
    
    #creates the B_hat matrix initialized to zeros
    B_hat = np.zeros((degree, degree))
    #iterates thorugh each degree
    for i in range(degree):

        #gets the degree for the b_d_M vector
        b_degree = (degree-i)

        #gets the b_d vector
        b_d = b_d_vector(degree=b_degree,
                         alpha=alpha)
        
        #gets the effective D
        D = np.eye(degree)
        #iterates through to get the complete D matrix
        for j in range(i):
            #gets the degree of the D matrix
            D_degree = degree - j
            #gets the next temp d matrix
            D_temp = D_d_M(d=D_degree,
                           M=0)
            #multiplies it into the D complete 
            D = D @ D_temp

        #now that we have the D matrix, we multiply it to the b vector
        b_vector_result = D @ b_d

        b_vector_result = b_vector_result.reshape((np.size(b_vector_result)))

        #puts it in the correct slot
        B_hat[:,i] = b_vector_result
    
    #gets  the inverse of the B_hat
    B_hat_inv = np.linalg.inv(B_hat)

    #returns both of the matrices
    return B_hat, B_hat_inv

#returns:
#1. initialized control points
#2. initialized status points
#creates matrix to initialize the control points
def initializeControlPoints(dimension: int,#the number of dimensions the spline exists within
                            degree: int, #the degree of the polynomials used in the Basis Spline
                            M: int): #the number of degrees of interest in the spline
    
    #sets the number of control points
    numControlPoints = M + degree
    
    #creates the vector to store the controlPoints
    controlPoints = np.full((dimension, numControlPoints), 0.0, dtype=float)

    #creates the vector of the metadata about the control points 
    # to determine whether a particular set of control points has been set by position and derivative conditions
    controlPointSetStatus = np.full((1, numControlPoints), False, dtype=bool)

    #returns the two vectors
    return controlPoints, controlPointSetStatus




#creates the function to perform the concatenation of B all together
def B_Concatenation(B_matrices: List[np.ndarray]) -> np.ndarray:

    #concatenates together the B matrices
    concatenatedMatrix = np.concatenate(B_matrices, axis=1)

    return concatenatedMatrix



#returns:
#1. the U1 matrix (corresponding to the linearly indepentent vectors)
#2. the U2 Matrix (corresponding to those in the null space)
#3. Sigma, the 2d Matrix of singular values
#4. V_T: the transpose V matrix (right eigenvectors)
#function to get the SVD of the full B matrix, which is the concatenation of the Tall
#B_M_d matrices
def B_SVD(B_complete: np.ndarray):#inputs the FULL B matrix. Not the B hat matrix mind you

    #calls the linear algebra svd

    U, Sigma_values, V_T = np.linalg.svd(B_complete)

    #gets the rank of B
    rank = np.linalg.matrix_rank(B_complete)
    #breaks it up into U1 and u2
    U1 = U[:,:rank]
    U2 = U[:,rank:]

    #creates the Sigma matrix
    Sigma = np.diag(Sigma_values)

    return U1, U2, Sigma, V_T


#function which creates the Q_M matrix
def Q_M(degree: int, #the degree of the bspline in question
        l: int, #the degree of the derivative we are working with
        M: int): #the number of intervals of interest in question
    
    apple = 0


#defines the function to create the S matrix
def S_M_d_l(M: int, #the number of intervals of interest in the particular spline
            degree: int, #the degree of the polynomials
            l: int): #the degree of the derivatives
    #starts out with the initial I matrix
    S_out = np.eye(M+degree)
    #iterates through d minus l
    for i in range(l+1):

        #gets the current D matrix using the function
        D_current = D_d_M(d=(degree - i), M=M)

        ##multiplies it into the S out matrix
        S_out = S_out @ D_current

    #returns the S out
    return S_out


#creates the function to get the integral of the b b transpose matrix
def integrate_b_bT(degree: int, #the degree of the bspline
                   l: int, #the degree of the derivative
                   M: int): #the number of intervals of interest
    #gets the length of the matrix
    matrixLength = M + degree - l
    #creates the zero matrix
    MatrixOutput = np.zeros((matrixLength, matrixLength))
    #iterates through the rows and columns
    for i in range(matrixLength):
        for j in range(matrixLength):
            #obtains the integration result
            integrationResult = integrate.quad(func=getb_bT_value, 
                                               a=0,
                                               b=M,
                                               args=(degree,1.0,l,M,i,j))
            #saves the integration result
            MatrixOutput[i,j] = integrationResult[0]

    
    #returns the matrix
    return MatrixOutput


#creates the function to calculate the W Matrix given M, d, l


#creates the function to get a specified value from the Matrix
def getb_bT_value(time: float,
                  degree: int,
                  alpha: float,
                  l: int,
                  M: int,
                  row: int,
                  col: int):
    
    #gets the degree of the b vector
    vectorDegree = degree - l
    #gets the vector from that
    vector = b_d_M_t_vector(time=time,
                            degree=vectorDegree,
                            alpha=alpha,
                            M=M)
    
    #gets the matrix
    A_matrix = vector @ np.transpose(vector)

    #gets the output value
    outputValue = A_matrix[row,col]

    #returns the output value
    return outputValue



#creates function to obtain the individual W matrix
def get_W_d_l_M(degree: int,
                 l: int,
                 M: int):
    
    #gets the S matrix
    S = S_M_d_l(M=M,degree=degree, l=(l-1))

    #calls the function to get the integrated matrix
    integrated_b_bT = integrate_b_bT(degree=degree, l=l, M=M)

    #gets the output Matrix
    W = S @ integrated_b_bT @ np.transpose(S)

    #returns the W matrix
    return W



#creates the function to get the complete W matrix
def get_W_d_M_rho(degree: int, #the degree of the polynomial
                  L: int, #the highest degree of the derivative (the highest it can be is d-1)
                  M: int,
                  rho: np.ndarray):
    
    #sets the number of iterations
    numIter = L + 1
    
    #asserts whether the input rho is of shape (L+1)x(1)
    assert rho.shape == (numIter,1), "position must have shape ({numIter},1}"


    #creates the zero matrix to begin the summation

    W_d_M_rho = np.zeros((M+degree, M+degree))
    #iterates through and obtains the corresponding W_d_l_M
    for l in range(numIter):
        
        #gets the temp W matrix
        temp_W = get_W_d_l_M(degree=degree, l=l, M=M)

        #adds the scales temp_w
        W_d_M_rho = W_d_M_rho + rho.item(l)*temp_W
        
        tornado = 0

    #returns the matrix
    return W_d_M_rho


#function to get a partition of the W complete matrix
def get_W_partitioned(degree: int,
                      M: int,
                      L: int,
                      rho: np.ndarray):
    
    #calls the function to get the W_d_M_rho matrix
    W_M = get_W_d_M_rho(degree=degree,
                        M=M,
                        L=L,
                        rho=rho)
    
    #gets the sections of it, and then 
    #gets the top row
    W_11 = W_M[0:degree, 0:degree]
    W_12 = W_M[0:degree, degree:M]
    W_13 = W_M[0:degree, M:(M+degree)]
    #gets the second row
    W_22 = W_M[degree:M, degree:M]
    W_23 = W_M[degree:M, M:(M+degree)]
    #gets the third row
    W_33 = W_M[M:(M+degree),M:(M+degree)]


    #gets the length of each partition
    length_1 = degree
    length_2 = M - degree
    length_3 = degree

    #reshapes all of them in the case that we get the whole shrinkage from 2d to 1d problem.
    #just in case.
    W_11 = getReshapedMatrix(M_in=W_11,
                             numRows=length_1,
                             numCols=length_1)
    
    W_12 = getReshapedMatrix(M_in=W_12,
                             numRows=length_1,
                             numCols=length_2)
    
    W_13 = getReshapedMatrix(M_in=W_13,
                             numRows=length_1,
                             numCols=length_3)
    
    W_22 = getReshapedMatrix(M_in=W_22,
                             numRows=length_2,
                             numCols=length_2)
    
    W_23 = getReshapedMatrix(M_in=W_23,
                             numRows=length_2,
                             numCols=length_3)
    
    W_33 = getReshapedMatrix(M_in=W_33,
                             numRows=length_3,
                             numCols=length_3)



    #gets the transpose components
    W_21 = np.transpose(W_12)
    W_31 = np.transpose(W_13)
    W_32 = np.transpose(W_23)

    #creates the top row
    row_1 = [W_11, W_12, W_13]
    #gets the second row
    row_2 = [W_21, W_22, W_23]
    #gets the third row
    row_3 = [W_31, W_32, W_33]

    #puts them together into a whole matrix
    W_parted = [row_1, row_2, row_3]

    #returns it 
    return W_parted





#defines a function to get reshape a matrix based on the length and the width
def getReshapedMatrix(M_in: np.ndarray, #the A matrix to be reshaped
                      numRows: int, #the number of rows of the matrix. the m in an mxn array
                      numCols: int): #the number of columns of the matrix. the n in and mxn array
    M_out = M_in.reshape((numRows, numCols))
    #returns the M_out
    return M_out


