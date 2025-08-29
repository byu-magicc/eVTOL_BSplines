#this file re implements the Path Generator's stuff, but this time rewritten to make more sense
#and to be a greater deal of optimality.


import os
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
import time




class SFC_PathGenerator:


    #creates the initialization function.
    def __init__(self,
                 num_intervals_free_space: int = 5):
        #sets the dimension to two for now
        self._dimension = 2
        #sets the number of intervals free space
        self.num_intervals_free_space = num_intervals_free_space
        #I'm only dealing with 3rd degree bsplines (called order in his stuff. I know - it's annoying)
        self._order = 3
        #gets the Minvo conversion matrix
        self._Minvo_matrix = get_M_matrix(order=self._order)
        #gets the objective functions from the c files
        self._objective_func_obj = ObjectiveFunctions(self._dimension)
        self._curvature_const_obj = CurvatureConstraints(self._dimension)
        self._waypoint_const_obj = WaypointConstraints(self._dimension)
        self._obstacle_cons_obj = ObstacleConstraints(self._dimension)



    #defines the function to generate a path through a list of SFCs with direct control point method.
    #that is, we directly modify and optimize the control points and see how that affects things
    def generatePath__directCtrlPts(self,
                                    numIntervalsOfInterestPerCorridor: int,
                                    initialControlPoints: np.ndarray,
                                    sfc_data: SFC_Data = None,
                                    objective_function_type: str = "minimal_velocity_path"):
        
        #gets the number of points



    #defines the function to generate the a path by modifying the start and end positions and the theta for velocity



#defines the function to get the number of control points from an existing array
def getNumCtrPts_array(controlPoints: np.ndarray):

    #gets the shape of the control points array
    controlPoints_shape = np.shape(controlPoints)

    #case it is along the zero axis
    if controlPoints_shape[0] > controlPoints_shape[1]:
        numControlPoints = controlPoints_shape[0]
    else:
        numControlPoints = controlPoints_shape[1]

    #returns the number of control points
    return numControlPoints


#gets the central control points from the complete initial control points list
def getPartitionedControlPoints_initial(initialControlPoints: np.ndarray,
                                        degree: int):
    
    #checks the shape and works with that shape
    pointsShape = np.shape(initialControlPoints)

    #case tall matrix
    if pointsShape[0] > pointsShape[1]:
        numControlPoints = pointsShape[0]
        numIntervalsOfInterest = numControlPoints - degree
        startControlPoints_initial = initialControlPoints[0:degree,:]
        centralControlPoints_initial = initialControlPoints[degree:numIntervalsOfInterest,:]
        endControlPoints_initial = initialControlPoints[numIntervalsOfInterest:,:]
    #case fat matrix
    else:
        numControlPoints = pointsShape[1]
        numIntervalsOfInterest = numControlPoints - degree
        startControlPoints_initial = initialControlPoints[:,0:degree]
        centralControlPoints_initial = initialControlPoints[:,degree:numIntervalsOfInterest]
        endControlPoints_initial = initialControlPoints[:,numIntervalsOfInterest:]

    #returns the initial central control points
    return startControlPoints_initial, centralControlPoints_initial, endControlPoints_initial

#function to reconstruct flattened control points, with the start and end sections
def reconstructFlattenedControlPoints(startControlPoints: np.ndarray,
                                      flattenedCenterControlPoints: np.ndarray,
                                      endControlPoints: np.ndarray,
                                      numCenterControlPoints: int,
                                      dimension: int):
    
    #reshapes the flattened center control points into the correct shape
    centerControlPoints = np.reshape(flattenedCenterControlPoints, (dimension, numCenterControlPoints))

    #concatenates them together to create the full control point array
    completeControlPoints = np.concatenate((startControlPoints, centerControlPoints, endControlPoints), axis=1)

    #returns  the control points
    return completeControlPoints