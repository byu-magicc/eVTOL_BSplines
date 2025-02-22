#This file implements the wrapper for the path generator submodule

#master class

import os, sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from bsplinegenerator.bsplines import BsplineEvaluation

#inserts the parent directory
sys.path.insert(0,os.fspath(Path(__file__).parents[2]))


#inserts the path needed for the path generator submodule
temp1 = os.fspath(Path(__file__).parents[1])
temp2 = os.path.abspath(os.path.join(temp1, 'src/submodules/path_generator'))
sys.path.insert(0,temp2)
temp3 = sys.path

tempPath = sys.path

#imports everything from the matrix evaluation 
from src.submodules.path_generator.path_generation.matrix_evaluation import *
from src.submodules.path_generator.path_generation.obstacle import *
from src.submodules.path_generator.path_generation.path_generator import *
from src.submodules.path_generator.path_generation.path_plotter import *
from src.submodules.path_generator.path_generation.safe_flight_corridor import *
from src.submodules.path_generator.path_generation.waypoint_data import *

#now that we have everything imported, how are we going to put it together?
#right now, it's all super duper clunky and fails the acid tests. 
#There has to be a solution. I just shotgunned imported everything.



#path_generator_simplified inherits from the original PathGenerator function
#in an attempt to simplify things down.
class path_generator_simplified(PathGenerator):


    #creates the init function
    #arguments:
    #1. dimension: int - the dimensionality we are working with
    def __init__(self,
                 dimension: int = 2):
        
        #calls the super init function for the original class
        super().__init__(dimension=dimension)


    
    #defines the waypoint path generator, which generates a path based solely on waypoints,
    #with no other constraints
    #Arguments:
    #1. waypoints: the given waypoints for the trajectory
    #2. max_curvature: the maximum allowable curvature (tightness of turns) for the thing
    #3. degree: the degree of the polynomial in question to be used
    #4. num_points_per_interval: the number of points output for the output spline data array 
    #5. max_incline: the maximum allowable incline, which we will set to none here I think
    def waypointPathGenerator(self,
                              waypoints: WaypointData,
                              max_curvature: float,\
                              degree: int,
                              num_points_per_interval: int,
                              max_incline: float = None)-> tuple[np.ndarray, np.ndarray, np.ndarray]:
        

        #gets the control points and the status using the generate path function
        controlPoints, status = self.generate_path(waypoint_data=waypoints,
                                                   max_curvature=max_curvature,
                                                   max_incline=max_incline,
                                                   sfc_data=None,
                                                   obstacles=None)
        
        #creates the bspline object from the control points
        bspline_object = BsplineEvaluation(control_points=controlPoints,
                                           order=degree,
                                           start_time=0.0,
                                           scale_factor=1,
                                           clamped=False)
        
        #creates the spline and time vector data
        spline_data, time_data = bspline_object.get_spline_data(num_data_points_per_interval=num_points_per_interval)
        
        #returns the spline_data, control points, and the time data
        return spline_data, controlPoints, time_data

#creates the wrapper 