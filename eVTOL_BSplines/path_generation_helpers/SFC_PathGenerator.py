#this is the file which implements the ability to generate 
#bspline paths through safe flight corridors. But this will
#be implemented differently from David Christensen's work, because
#it will be done by 

import numpy as np

#the importation of the Safe Flight Corridor, and paths of Safe Flight Corridors
from submodules.path_generator.path_generation.safe_flight_corridor import SFC, SFC_Data
#the importation of the waypoint, which refers to the position, acceleration, velocity, etc at a given time
#and then WaypointData refers to a list of those, so like beginning and end positions with intermediate positions
from submodules.path_generator.path_generation.waypoint_data import Waypoint, WaypointData
#imports the function to get the matrix M which converts from standard control points to minvo control points
from submodules.path_generator.path_generation.matrix_evaluation import get_M_matrix, evaluate_point_on_interval




import time


#creates the PathGenerator Class. This is really frustrating because I don't know exactly what's
#supposed to be changed and what's the best way to make this work, and what can I modify?
#and what should I not modify?
class PathGenerator:


    #creates the initialization function
    def __init__(self,
                 dimension: int,
                 num_intervals_free_space: int = 5,
                 degree: int = 3):
        
        #saves the dimension of this problem
        self.dimension = dimension
        #saves the degree of the polynomial B-Spline
        self.degree = degree
        #saves the M matrix, which is the Minvo matrix
        self.M_d = get_M_matrix(order=degree)