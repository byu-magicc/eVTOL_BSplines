#This file is created to ease the process of plotting out the paths
#before this, all the paths were super disjointed and impossible to work with
#the goal here is to fix that at least partially


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


#imports the function from the other path plotter, which all it does is keep the axes constant
from eVTOL_BSplines.submodules.path_generator.path_generation.path_plotter import set_axes_equal





#creates the pathplotter class
class pathPlotter:

    #creates the initialization function
    def __init__(self):

        