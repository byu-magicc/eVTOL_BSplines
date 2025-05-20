#%%
import os, sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from bsplinegenerator.bsplines import BsplineEvaluation
from IPython.display import display
import sympy as sp
import time

from IPython.display import display

import scipy.integrate as integrate

import matplotlib.pyplot as plt

sys.path.insert(0,os.fspath(Path(__file__).parents[1]))
tempPath = sys.path


from eVTOL_BSplines.path_generation_helpers.control_points_matrix import controlPointsGenerator

degree = 3
M=50
dimension = 2

rho = np.array([[1.0,1.0,1.0]])

ctrlPntsGen = controlPointsGenerator(dimension=dimension,
                                     degree=degree,
                                     M=M)

#creates the initial and final conditions
Start = np.array([[0.0, 0.0, 0.0],
                  [1.0, 1.0, 1.0]])


End = np.array([[10.0, 10.0, 1.0],
                [5.0, 0.0, 0.0]])

print("Start Conditions")
display(sp.Matrix(Start))
print("End Conditions")
display(sp.Matrix(End))


start_time = time.time()
#gets the control points
ctrlPnts = ctrlPntsGen.generateControlPoints(Start=Start,
                                             End=End,
                                             rho=rho)
end_time = time.time()

function_time_delay = end_time - start_time
print(f"Control Points Time Delay: {function_time_delay} seconds")

#evaluates with the control points
evaluation = BsplineEvaluation(control_points=ctrlPnts,
                               order=degree,
                               start_time=0.0,
                               scale_factor=1.0,
                               clamped=False)

#gets the spline and time data
splineData, timeData = evaluation.get_spline_data(num_data_points_per_interval=20)

plt.figure(0)
plt.title("Spline Evaluation")
plt.plot(splineData[0,:], splineData[1,:], linewidth=3.0)
#plots the first arrow
plt.arrow(x=Start[0,0],
          y=Start[1,0],
          dx=Start[0,1],
          dy=Start[1,1],
          color='green')
plt.arrow(x=End[0,0],
          y=End[1,0],
          dx=End[0,1],
          dy=End[1,1],
          color='red')
plt.show()



potato = 0
# %%
