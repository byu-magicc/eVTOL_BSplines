#this file is for helping out with the matrix representation. 
#I wish that some of the functionality had been included in the path generator files. 
#alas, they were not, but that's okay for now.

import os, sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from bsplinegenerator.bsplines import BsplineEvaluation

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
                                 scale_factor: float, #the factor of scaling for the knot vector
                                 start_time: float): #the start time for the zeroth knot point
    
    #gets the number of knot points
    num_knot_points = M + 2*degree + 1

    #gets the knot points, which is in the range of zero to the number of knot points
    #then it is shifted back by the degree
    #and then it is scaled by the scaling factor
    #finally, the start time defines the initial start time shift
    knot_points = (np.arange(num_knot_points) - degree)*scale_factor + start_time


    #returns the knot_points
    return knot_points


#creates a uniform bspline wrapper for the cox de boor function
def uniform_cox_de_Boor_basis_function(time: float, #the time at which to evaluate the function
                                       degree: int, #the degree of the function to evaluate
                                       knot_points: np.ndarray): #the knot vector

    
    #gets the largest index s where knotPoints[s] < t
    #TODO make sure the minus one is actually working correctly
    s = np.searchsorted(knot_points, time) - 1

    #gets the length of the knot points
    numKnotPoints = np.size(knot_points)
    #gets the index for the end time
    end_time_index = int(numKnotPoints - (degree + 1))
    #gets the end time
    end_time = knot_points.item(end_time_index)

    #gets the particular table
    table = cox_de_boor_table_basis_function(time=time,
                                             i=s,
                                             order=degree,
                                             knot_points=knot_points,
                                             end_time=end_time,
                                             clamped=False)
    
    #returns the table
    return table

#defines the function to get the uniform cox de Boor basis table
def uniform_cox_de_boor_basis_function_table(time: float, #the time at which to evalues the functions
                                             degree: int, #the degree of the basis functions
                                             knot_points: np.ndarray): #the knot points for the vector
    
    #gets the largest index s where knotPoints[s] < t
    #TODO make sure the minus one is actually working correctly
    s = np.searchsorted(knot_points, time) - 1

    #gets the length of the knot points
    numKnotPoints = np.size(knot_points)
    #gets the index for the end time
    end_time_index = int(numKnotPoints - (degree + 1))
    #gets the end time
    end_time = knot_points.item(end_time_index)
    #creates the table 
    table = __cox_de_boor_table_basis_function_whole_table(time=time,
                                                           i=s,
                                                           order=degree,
                                                           knot_points=knot_points,
                                                           end_time=end_time,
                                                           clamped=False)
    
    #returns the table
    return table



#function which returns the value of a basis function at a particular time 
#gets the vector b given degree d, length M, and particular evaluation time t
def b_d_M_t(t: float, #particular time of evaluation
            d: int, #degree of the polynomial to evaluate
            M: int): #related to length of evaluation
    
    #the length of the vector is d + M
    vectorLength = d + M

    #creates the vector of zeros
    b_d_M = np.zeros((vectorLength, 1))

    for i in range(vectorLength):
        b_d_temp = 0
        #gets  the evaluation of the particular bspline
        b_d_M[i,0] = 0




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


#function that gets the B_d_M matrix
def get_B_d_M(b_d_M: np.ndarray, #the b vector for the evaluation of the basis functions
              d: int, # the degree of the polynomial
              M: int): # the length of the interval of interest

    #creates the initial matrix

    #there are d+1 columns in the matrix
    for i in range(d + 1):
        potato = 0

#creates the functions for degrees 1 through 5