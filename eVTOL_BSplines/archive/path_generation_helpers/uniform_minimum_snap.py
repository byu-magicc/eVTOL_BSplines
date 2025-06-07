



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


from eVTOL_BSplines.path_generation_helpers.conditions_helpers import conditions, conditionsList



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

    #creates the initialization function
    def __init__(self,
                 listOfConditions: conditionsList): #inputs the conditions list using the 
        
        #saves the list of conditions
        
