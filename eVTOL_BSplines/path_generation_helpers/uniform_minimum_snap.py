



import os, sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from bsplinegenerator.bsplines import BsplineEvaluation
from scipy.optimize import minimize, Bounds, LinearConstraint, NonlinearConstraint
import math
from enum import Enum


sys.path.insert(0,os.fspath(Path(__file__).parents[2]))

temp1 = os.fspath(Path(__file__).parents[0])
temp2 = os.path.abspath(os.path.join(temp1, 'submodules/path_generator'))
sys.path.insert(0,temp2)
tempPath = sys.path


from eVTOL_BSplines.path_generation_helpers.conditions_helpers import conditions



#creates the objective type for the trajectory generator
class ObjectiveType(Enum):
    MINIMIZE_TIME = 1
    MINIMIZE_DISTANCE = 2
    MINIMIZE_TIME_AND_DISTANCE = 3
    MINIMIZE_VELOCITY = 4
    MINIMIZE_ACCELERATION = 5
    MINIMIZE_JERK = 6
    MINIMIZE_SNAP = 7


#creates the constrained trajectory generator function
class UniformConstrainedTrajectory:

    #creeates the initialization function
    def __init__(self,
                 objectiveType: ObjectiveType,
                 dimension: int = 2,
                 max_interval_distance = 1,
                 control_point_bounds = [-100,100],
                 order: int = 5):
        
        #saves the information
        self._objectiveType = objectiveType
        self._max_interval_distance = max_interval_distance
        self._control_point_bounds = control_point_bounds
        self._dimension = dimension
        self._order = order
        self._max_curvature = 10000


    #creates the generator trajectory function
    def generate_trajectory(self, waypoints, max_curvature=np.inf):
        if self._objectiveType == ObjectiveType.MINIMIZE_SNAP:
            objectiveFunction = self.__minimize_snap_objective_function
        self._max_curvature = max_curvature
        initial_control_points = self.__create_interpolated_points(waypoints)
        number_of_control_points = np.shape(initial_control_points)[1]
        print("number_control_points: " , number_of_control_points )
        initial_scale_factor = 1
        optimization_variables = np.concatenate((initial_control_points.flatten(),[initial_scale_factor]))
        optimization_variable_lower_bound = optimization_variables*0 + self._control_point_bounds[0]
        optimization_variable_upper_bound = optimization_variables*0 + self._control_point_bounds[1]
        optimization_variable_lower_bound[-1] = 0.0001
        optimization_variable_bounds = Bounds(lb=optimization_variable_lower_bound, ub = optimization_variable_upper_bound)
        waypoint_constraint = self.__create_waypoint_constraint(waypoints, number_of_control_points)
        result = minimize(
            objectiveFunction,
            x0=optimization_variables,
            method='SLSQP',
            bounds = optimization_variable_bounds,
            constraints=(waypoint_constraint))
        optimized_scale_factor = result.x[-1]
        control_points_optimized = result.x[0:number_of_control_points*self._dimension].reshape(self._dimension,number_of_control_points)
        return control_points_optimized, optimized_scale_factor




    #creates the function to obtain the minimum snap trajectory
    def minimumSnapControlPoints(self, )