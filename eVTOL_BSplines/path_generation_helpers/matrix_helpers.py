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
def get_D_d_M(d: int,  #the degree of the spline
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
def B_d_M_t_matrix(time: float, #time at which to evaluate the matrix
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
            D_temp = get_D_d_M(d=D_degree,
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
    return B_d_M, B_hat_d_M


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
            D_temp = get_D_d_M(d=D_degree,
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




