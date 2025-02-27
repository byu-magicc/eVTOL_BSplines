import os, sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from bsplinegenerator.bsplines import BsplineEvaluation

sys.path.insert(0,os.fspath(Path(__file__).parents[1]))

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



#creates the waypoint path generator class
class waypointPathGenerator(PathGenerator):

    #creates the init function
    #arguments:
    #1. dimension: int - the dimensionality we are working with
    #2. the waypoints in the form of the waypoint data class
    def __init__(self,
                 waypoints: WaypointData,
                 max_curvature: float,
                 degree: int,
                 num_points_per_interval: int,
                 max_incline: float = None,
                 dimension: int = 2,
                 ):
        
        #calls the super init function for the original class
        super().__init__(dimension=dimension)

        self.num_data_points_per_interval = num_points_per_interval

        #calls the function to generate the control points based on the constraints
        self.controlPoints, status = self.generate_path(waypoint_data=waypoints,
                                                   max_curvature=max_curvature,
                                                   max_incline=max_incline,
                                                   sfc_data=None,
                                                   obstacles=None)
        
        #creates the bspline object based on this data
        self.bspline_object = BsplineEvaluation(control_points=self.controlPoints,
                                                order=degree,
                                                start_time=0.0,
                                                scale_factor=1,
                                                clamped=False)
        

    #gets the control points
    def getControlPoints(self)->np.ndarray:
        return self.controlPoints



    
    #creates the function to get the positional data
    def getPosData(self)->tuple[np.ndarray, np.ndarray]:
        #calls the get positional and time data functions and returns them now
        pos_data, pos_time_data = self.bspline_object.get_spline_data(num_data_points_per_interval=self.num_data_points_per_interval)
        #returns that data in the function
        return pos_data, pos_time_data
    

    #creates the function to get the velocity data
    def getVelData(self)->tuple[np.ndarray, np.ndarray]:

        #calls the funciton to get the velocity and time data functions
        vel_data, vel_time_data = self.bspline_object.get_spline_derivative_data(num_data_points_per_interval=self.num_data_points_per_interval,
                                                                                 rth_derivative=1)
        
        #returns the data
        return vel_data, vel_time_data

    #creates the function to get the accel data
    def getAccelData(self)->tuple[np.ndarray, np.ndarray]:

        #calls the funciton to get the velocity and time data functions
        accel_data, accel_time_data = self.bspline_object.get_spline_derivative_data(num_data_points_per_interval=self.num_data_points_per_interval,
                                                                                     rth_derivative=2)
        
        #returns the data
        return accel_data, accel_time_data