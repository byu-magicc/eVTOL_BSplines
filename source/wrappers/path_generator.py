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
temp2 = os.path.abspath(os.path.join(temp1, 'submodules/path_generator'))
sys.path.insert(0,temp2)
temp3 = sys.path

tempPath = sys.path

#imports everything from the matrix evaluation 
from source.submodules.path_generator.path_generation.matrix_evaluation import *
from source.submodules.path_generator.path_generation.obstacle import *
from source.submodules.path_generator.path_generation.path_generator import *
from source.submodules.path_generator.path_generation.path_plotter import *
from source.submodules.path_generator.path_generation.safe_flight_corridor import *
from source.submodules.path_generator.path_generation.waypoint_data import *

#now that we have everything imported, how are we going to put it together?
#right now, it's all super duper clunky and fails the acid tests. 
#There has to be a solution. I just shotgunned imported everything.



#creates the wrapper 