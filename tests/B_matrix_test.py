#I'm just using this file to test the B matrices generation


import os, sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from bsplinegenerator.bsplines import BsplineEvaluation
from IPython.display import display

sys.path.insert(0,os.fspath(Path(__file__).parents[1]))
tempPath = sys.path



#imports the bspline
from eVTOL_BSplines.path_generation_helpers.matrix_helpers import B_M_matrix


M = 5
d = 4


#gets the B matrix for time = M, and time = 0
B_0 = B_M_matrix(time=0.0,
                 degree=d,
                 alpha=1.0,
                 M=M)


B_M = B_M_matrix(time=M,
                 degree=d,
                 alpha=1.0,
                 M=M)




potato = 0