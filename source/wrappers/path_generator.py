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

from source.submodules.path_generator.path_generation.path_generator import PathGenerator
from source.submodules.path_generator.path_generation.safe_flight_corridor import SFC, SFC_Data, plot_sfcs, get2DRotationAndTranslationFromPoints
from source.submodules.path_generator.path_generation.obstacle import Obstacle, plot_2D_obstacles
from source.submodules.path_generator.path_generation.waypoint_data import Waypoint, WaypointData
from source.submodules.path_generator.path_generation.path_plotter import set_axes_equal
import time

