import os, sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from bsplinegenerator.bsplines import BsplineEvaluation
from IPython.display import display
import sympy as sp

from IPython.display import display

import scipy.integrate as integrate

sys.path.insert(0,os.fspath(Path(__file__).parents[1]))
tempPath = sys.path


degree = 5
M=50
rho = np.array([[1.0,1.0,1.0,1.0,1.0]])

from eVTOL_BSplines.path_generation_helpers.control_points_generator import PathGenerator



#creates the start and end conditions
S = np.array([[1,1,1,1,1],
              [1,1,1,1,1]])

E = np.array([[1,1,1,1,1],
              [1,1,1,1,1]])


pathGen = PathGenerator(dimension=2,
                        M=M,
                        degree=degree)

#gets the control points
pathGen.getControlPoints(S=S,
                         E=E,
                         rho=rho)



