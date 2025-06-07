#this file creates the implementation to create a spline
#using the matrix method of obtaining control points
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

from eVTOL_BSplines.path_generation_helpers.matrix_helpers import *
from eVTOL_BSplines.path_generation_helpers.matrix_generators_efficient import *

from eVTOL_BSplines.submodules.path_generator.path_generation.waypoint_data import *


class splineGeneratorMatrix:

    #creates the initialization function
    def __init__(self,
                 waypoints: WaypointData,
                 max_curvature: float,
                 degree: int,
                 num_points_per_interval: int,
                 max_incline: float,
                 dimension: int = 2):

        #calls the super init function on the original class for Pathgenerator
        super().__init__(dimension=dimension)
        
        self.waypointsData = waypoints

        self.num_data_points_per_interval = num_points_per_interval







        
