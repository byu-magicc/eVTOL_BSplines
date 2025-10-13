#this file re implements the Path Generator's stuff, but this time rewritten to make more sense
#and to be a greater deal of optimality.

import os, sys
import numpy as np
from scipy.optimize import minimize, Bounds, LinearConstraint, NonlinearConstraint, Bounds
from path_generation.matrix_evaluation import get_M_matrix, evaluate_point_on_interval
from PathObjectivesAndConstraints.python_wrappers.objective_functions import ObjectiveFunctions
from PathObjectivesAndConstraints.python_wrappers.curvature_constraints import CurvatureConstraints
from PathObjectivesAndConstraints.python_wrappers.obstacle_constraints import ObstacleConstraints
from PathObjectivesAndConstraints.python_wrappers.incline_constraints import InclineConstraints
from PathObjectivesAndConstraints.python_wrappers.waypoint_constraints import WaypointConstraints
from bsplinegenerator.bspline_to_minvo import get_composite_bspline_to_minvo_conversion_matrix
from path_generation.safe_flight_corridor import SFC_Data, SFC
from path_generation.obstacle import Obstacle
from path_generation.waypoint_data import Waypoint, WaypointData
from eVTOL_BSplines.message_types.msg_control_points import MSG_Control_Points
from eVTOL_BSplines.message_types.msg_annulus_convex_hull import Msg_Annulus_Convex_Hull
import time
import cProfile
import pstats
from copy import deepcopy
from bsplinegenerator.bspline_to_minvo import convert_list_to_minvo_control_points
from path_generation.sfc_path_generator_helpers import *

import cvxpy as cp


class SFC_PathGenerator:

    def __init__(self,
                 dimension: int = 3,
                 num_intervals_free_space: int = 5,
                 degree: int = 3,
                 M: int = 10):

        self.dimension = dimension
        self.num_intervals_free_space = num_intervals_free_space
        self.degree = degree
        self.M = M

        self._Minvo_matrix = get_M_matrix(order=self.degree)

    #defines the function to generate the path using the control points
    def generatePath(self,
                     startControlPoints: np.ndarray,
                     endControlPoints: np.ndarray,
                     sfc_data: SFC_Data,
                     numPointsPerUnit: float):
        
        #based on the SFC_Data, we get the list of number of control points designated for each 
        numCntPts_list = getNumCntPts_list(sfc_data=sfc_data,
                                           numPointsPerUnit=numPointsPerUnit)

        #gets the total number of contorl points from the list of number of control points
        numCntPts_total = getNumCntPts(numCntPts_list=numCntPts_list,
                                       degree=self.degree)

        #instantiates the control points variable.
        controlPoints_cpVar = cp.Variable((self.dimension, numCntPts_total))


        #calls the get overlapping control points constraints function
        controlPoints_constraints = \
        self.get_overlapping_control_points_constraints(controlPoints_var=controlPoints_cpVar,
                                                        sfc_data=sfc_data,
                                                        numCntPts_list=numCntPts_list,
                                                        startControlPoints=startControlPoints,
                                                        endControlPoints=endControlPoints)
        

        #obtains the objective function
        objective_function = self.minimum_distance_objective(cpVar_cntPts=controlPoints_cpVar)



        #############################################################################
        #actual solving problem section
        problem = cp.Problem(objective=objective_function,
                             constraints=controlPoints_constraints)
        
        #calls the problem solver function
        problem.solve(solver=cp.CLARABEL)

        outputControlPoints = controlPoints_cpVar.value

        #returns the output contorl points
        return outputControlPoints



    #defines the function to create overlapping control point constraints
    def get_overlapping_control_points_constraints(self,
                                                   controlPoints_var: cp.Variable,
                                                   sfc_data: SFC_Data,
                                                   numCntPts_list: list[int],
                                                   startControlPoints: np.ndarray,
                                                   endControlPoints: np.ndarray):
        
        #gets the number of dimensions
        controlPoints_shape = controlPoints_var._shape
        #and the number of dimensions
        numDimensions = controlPoints_shape[0]

        totalStartTime = time.time()
        #gets the sfc list
        sfc_list = sfc_data.get_sfc_list()

        #creates the control poitns index
        controlPoints_index = 0 

        #creates the list of conditions
        controlPoints_constraints = []

        #creates the list of lists for the A matrices for each control point
        A_matrices_SFC_list = []
        #does the same for the b vectors
        b_vectors_sfc_list = []

        A_matrices_complete = []

        b_vectors_complete = []

        #General SFC constraints section
        #iterates over each safe flight corridor in the sfc list
        for i, sfc in enumerate(sfc_list):

            #gets the number of control points in the current corridor
            numCntPoints_inCorridor = numCntPts_list[i]

            #gets the incremental index
            incremental_index = numCntPoints_inCorridor - self.degree

            #gets the current control points partition
            controlPointsPartition \
                  = controlPoints_var[:,(controlPoints_index):(controlPoints_index+numCntPoints_inCorridor)]

            #gets the normals and vertices
            normals, vertices = sfc.getNormalsVertices()

            #gets the A and b matrices
            A, b = generate_A_b(normalVectors=normals,
                                vertices=vertices)
            

            for i in range(numCntPoints_inCorridor):
                #current index
                currentIndex = controlPoints_index + i

                #checks if the index currently exists in the A list
                if currentIndex < len(A_matrices_complete):
                    #then we can add to it
                    #gets the current A and b lists
                    A_matrices_complete[currentIndex].append(A)
                    b_vectors_complete[currentIndex].append(b)
                #otherwise we create it
                else:
                    A_matrices_complete.append([A])
                    b_vectors_complete.append([b])
            
            A_matrices_SFC_list.append(A)
            b_vectors_sfc_list.append(b)

            #appends the A and b matrices to the respective sections in the list
            
            #gets the temp inequality constraint
            inequalityConstraint_temp = [A @ controlPointsPartition <= b]

            #adds this to the constraints
            controlPoints_constraints += inequalityConstraint_temp

            #increments the control points index by the incremental amount
            controlPoints_index += incremental_index


        #creates the concatenated A and b matrices
        self.A_cat_list = []

        self.b_cat_list = []

        for A_list, b_list in zip(A_matrices_complete, b_vectors_complete):

            A = np.ndarray((0,numDimensions))
            b = np.ndarray((0,1))
            #iterates over the items here
            for A_temp, b_temp in zip(A_list, b_list):
                A = np.concatenate((A, A_temp), axis=0)
                b = np.concatenate((b, b_temp), axis=0)

            #adds this A and b to the full list
            self.A_cat_list.append(A)
            self.b_cat_list.append(b)


        #gets the varaibles corresponding to the start and end control points
        startControlPoints_cpVar = controlPoints_var[:, :self.degree]
        endControlPoints_cpVar = controlPoints_var[:, (-self.degree):]


        #creates the equality constraints for start and end conditions
        startEqualityConstraint = [startControlPoints_cpVar == startControlPoints]
        endEqualityConstraint = [endControlPoints_cpVar == endControlPoints]

        #adds these equality constraints to the main constraints list
        controlPoints_constraints += startEqualityConstraint
        controlPoints_constraints += endEqualityConstraint

        totalEndTime = time.time()
        totalTime = totalEndTime - totalStartTime

        #returns the constraints list
        return controlPoints_constraints


    #creates the function to define the objective function
    def minimum_distance_objective(self, cpVar_cntPts: cp.Variable):
        #gets the velocity control poitns
        velocityControlPoints_cp = cpVar_cntPts[:,0:-1] - cpVar_cntPts[:,1:]

        minimizeLength_objectiveFunction = cp.Minimize(cp.sum(cp.norm(velocityControlPoints_cp, axis=1)))

        #returns the objective
        return minimizeLength_objectiveFunction


#defines the function to get the absolute index of the control points, given the indexing list
def getAbsoluteIndex(numCntPts_list: list[int],
                     sfc_index: int,
                     controlPoint_index: int,
                     degree: int):
    
    tempSum = 0
    #iterates to get the previous sum
    for i in range(sfc_index):
        tempSum += (numCntPts_list[i] - degree)

    #adds the temp sum to the contorl point index
    tempSum += controlPoint_index

    return tempSum

#defines to get the SFC indices from the control point absolute index
def getSFCIndices_absoluteIndex(controlPoint_absoluteIndex: int,
                                sfcLengths_list: list[int],
                                degree: int):

    #gets the number of corridors
    numCorridors = len(sfcLengths_list)

    sum = 0

    for i, sfc_length in enumerate(sfcLengths_list):
        
        #case we are on the zeroth one
        if i == 0:

            #case we are in the start or center portion
            if sum <= controlPoint_absoluteIndex < (sum + sfc_length - degree):
                sfc_indices = [i]
                return sfc_indices
            #case we are in the overlap
            elif (sum + sfc_length - degree) <= controlPoint_absoluteIndex < (sum + sfc_length):
                #then we return two
                sfc_indices = [i, i+1]
                return sfc_indices
        
        #case we are in the middle ones
        elif 0 < i < (numCorridors - 1):

            #case we are in the beginning overlap section
            if sum <= controlPoint_absoluteIndex < (sum + degree):
                sfc_indices = [i-1, i]
                return sfc_indices
            elif (sum + degree) <= controlPoint_absoluteIndex < (sum + sfc_length - degree):
                sfc_indices = [i]
                return sfc_indices
            elif (sum + sfc_length - degree) <= controlPoint_absoluteIndex < (sum + sfc_length):
                #then we return two
                sfc_indices = [i, i+1]
                return sfc_indices
        elif i == (numCorridors - 1):
            #case we are in the beginning overlap section
            if sum <= controlPoint_absoluteIndex < (sum + degree):
                sfc_indices = [i-1, i]
                return sfc_indices
            elif (sum + degree) <= controlPoint_absoluteIndex < (sum + sfc_length - degree):
                sfc_indices = [i]
                return sfc_indices
        
        #then increments by the sfc length minus the degree
        sum += (sfc_length - degree)



#arguments:
#M: the number of control points in the SFC
#N: the number of control points used in the Annulus
#d: the degree of the bspline
def getLowerIndex(numCntPts_SFC: int,
                  numCntPts_Annulus: int,
                  d: int):

    bottomIndex = int(round(numCntPts_SFC - ((numCntPts_Annulus+d)/2)))

    return bottomIndex

#gets the upper index
def getUpperIndex(numCntPts_Annulus: int,
                  d: int):
    topIndex = int(round((numCntPts_Annulus + d)/2)) - 1

    return topIndex