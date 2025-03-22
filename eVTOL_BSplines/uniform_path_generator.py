#This file implements a b spline creation algorithm which makes a b spline, and which follows Dr. Beard's specifications for the design

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


#creates the uniform path generator
class uniformPathGenerator:

    #creates the initialization function
    def __init__(self,
                 dimension: int,
                 num_intervals_free_space: int = 5):

        #saves the dimension of the problem
        self.dimension = dimension
        #saves the number of intervals for free space
        self.num_intervals_free_space = num_intervals_free_space

    #creates the generate path function
    def generate_path(self, initial_conditions: )
