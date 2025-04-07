#%%

import os, sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from bsplinegenerator.bsplines import BsplineEvaluation
from IPython.display import display

from bsplinegenerator.bsplines import BsplineEvaluation
from bsplinegenerator.bspline_plotter import plot_bspline

sys.path.insert(0,os.fspath(Path(__file__).parents[1]))
tempPath = sys.path



#imports the bspline
from eVTOL_BSplines.path_generation_helpers.matrix_helpers import getCtrlPtswSVD

degree = 3
M = 15
alpha = 1.0


#creates the Start Conditions
S = np.array([[0, 0, 0, 0, 0],
              [0.1, 0.1, -9.81, 0, 0]])

#sets the end conditions
E = np.array([[100, 20, 0.1, 0, 0],
              [-10,   0,  0, 0, 0]])


#creates the rho
rho = np.array([[0.0],[0.0],[0.0],[0.0],[0.0]])

#gets the control points
C = getCtrlPtswSVD(S=S,
                   E=E,
                   degree=degree,
                   M=M,
                   alpha=alpha,
                   rho=rho)

display(C)


bspline = BsplineEvaluation(control_points=C,
                            order=degree,
                            start_time=0.0,
                            scale_factor=1.0,
                            clamped=False)

bspline.plot_spline(num_data_points_per_interval=10)



potato = 0

# %%
