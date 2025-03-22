#This file is created to ease the process of plotting out the paths
#before this, all the paths were super disjointed and impossible to work with
#the goal here is to fix that at least partially


import os, sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from bsplinegenerator.bsplines import BsplineEvaluation
from typing import Optional


sys.path.insert(0,os.fspath(Path(__file__).parents[1]))

temp1 = os.fspath(Path(__file__).parents[0])
temp2 = os.path.abspath(os.path.join(temp1, 'submodules/path_generator'))
sys.path.insert(0,temp2)
tempPath = sys.path


#imports the function from the other path plotter, which all it does is keep the axes constant
from eVTOL_BSplines.submodules.path_generator.path_generation.path_plotter import set_axes_equal
from eVTOL_BSplines.submodules.path_generator.path_generation.waypoint_data import Waypoint, WaypointData




#creates the function to plot a bspline, 
# and optionally to plot waypoints and control points
#arguments:
#1. bspline data: the positional data of the bspline,
#2. waypoints: the waypoints, which is of the class Waypoint data (optional)
#3. control_points: the location of the control points (optional)
def plotSpline_2d(bspline_data: np.ndarray,
                  waypoints: Optional[WaypointData] = None,
                  control_points: Optional[np.ndarray] = None):
    
    #gets the dimension of the bspline data

    ax = plt.axes()
    ax.plot(bspline_data[0,:], bspline_data[1,:])

    #checks if waypoints is not of None type
    if waypoints is not None:
        #gets the waypoints locations
        waypoint_locations = waypoints.get_waypoint_locations()
        #plots the waypoints otherwise
        ax.scatter(waypoint_locations[0,:], waypoint_locations[1,:], c='blue')

    #checks if the controlPoints is not of the none type
    if control_points is not None:
        #scatter plots the control points
        ax.scatter(control_points[0,:], control_points[1,:], c='orange')


    set_axes_equal(ax=ax, dimension=2)
    plt.show()


#function to plot the spline, the actual position of the craft, and the calculated error
#Arguments:
#1. bspline_data: a 2d array of the bspline plot data (usually the commanded positions)
#2. position_data: a 2d array of the actual aircraft position through the thing
#3. plot error: bool to set whether we are plotting the error, or the difference between the two
def plotSplinePositionerror_2d(bspline_data: np.ndarray,
                               position_data: np.ndarray,
                               plot_error: bool = False,
                               waypoints: Optional[WaypointData] = None):
    
    ax = plt.axes()
    ax.plot(bspline_data[0,:], bspline_data[1,:], c='blue', label='BSpline')
    ax.plot(position_data[0,:], position_data[1,:], c='orange', label='Position Data')

    #gets the error and plots it out as well, if plot error is true
    if plot_error:
        pos_error = bspline_data - position_data
        ax.plot(pos_error[0,:], pos_error[1,:], c='green', label='Position Error')


    #checks if waypoints is not of None type
    if waypoints is not None:
        #gets the waypoints locations
        waypoint_locations = waypoints.get_waypoint_locations()
        #plots the waypoints otherwise
        ax.scatter(waypoint_locations[0,:], waypoint_locations[1,:], c='blue')

    set_axes_equal(ax=ax, dimension=2)
    ax.legend()
    plt.show()
    

