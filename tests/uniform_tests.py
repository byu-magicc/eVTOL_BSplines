import os, sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from bsplinegenerator.bsplines import BsplineEvaluation

sys.path.insert(0,os.fspath(Path(__file__).parents[1]))
tempPath = sys.path

from eVTOL_BSplines.path_generation_helpers.matrix_helpers import uniform_cox_de_Boor_basis_function, uniform_knot_point_generator, uniform_cox_de_boor_basis_function_table
from eVTOL_BSplines.path_generation_helpers.matrix_helpers import uniform_basis_function_evaluation


#sets the number of intervals of interest
M = 9
#sets the degree of the spline
degree = 5
#sets the number of control points
Num_Ctrl_Pts = M + degree

#defines the initial acceleration to be gravity
accel_init = -9.81


#evaluates the uniform basis function at a time
evaluation = uniform_basis_function_evaluation(time=1.5,degree=2)
print(evaluation)

numTimeSteps = 400
#creates a time vector
times = np.linspace(start=-1.0, stop=6.0, num=numTimeSteps)
evaluations = []
for i in range(numTimeSteps):
    temp = uniform_basis_function_evaluation(time=times[i], degree=1)
    evaluations.append(temp)

evaluations = np.array(evaluations)

plt.figure(0)
plt.plot(evaluations)
plt.show()


#gets the knot point vector
knotPoints = uniform_knot_point_generator(M=M,
                                          degree = degree,
                                          alpha=1,
                                          start_time=0)


#evaluates the uniform de boor basis function at time = 0
value = uniform_cox_de_Boor_basis_function(time=0.0,
                                           degree=degree,
                                           knot_points=knotPoints)


#gets the whole table
table = uniform_cox_de_boor_basis_function_table(time=0.0,
                                                 degree=degree,
                                                 knot_points=knotPoints)


#creates the initial P conditions
P_0 = np.array([[0, 0, 0, 0, 0],
                [0,0,accel_init, 0, 0]])


#sets the horizontal x distance to travel for the transition
x_distance = 40.0
#sets the altitude
cruise_altitude = 10.0

#sets the cruise velocity
cruise_velocity = 25.0


#creates the final P conditions
P_M = np.array([[x_distance, cruise_velocity, 0, 0, 0],
                [-cruise_altitude, 0, 0, 0, 0]])




potato = 0