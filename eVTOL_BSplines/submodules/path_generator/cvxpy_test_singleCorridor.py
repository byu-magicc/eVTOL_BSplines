#this file is my first preliminary test to see how cvxpy actually works
import cvxpy as cp
import numpy as np
from numpy import cos, sin


from bsplinegenerator.bsplines import BsplineEvaluation

import matplotlib.pyplot as plt



#
A = np.array([[-1, 0],
              [0, 1],
              [1, 0],
              [0, -1]])



b = np.array([[-1],
              [6],
              [3],
              [-1]])



M = 10

d = 3

numControlPoints = M + d

#creates the initial and final control points
initialControlPoints = np.array([[2, 2, 2],
                                 [1.25, 1.5, 1.75]])

finalControlPoints = np.array([[2.25, 2.0, 1.75],
                               [5.25, 5.5, 5.75]])



#creates the cp variables
controlPoints_cvxpyVar = cp.Variable((2, numControlPoints))


controlPoints_constraints = []


boundaryConstraint = [A @ controlPoints_cvxpyVar <= b]


controlPoints_constraints += boundaryConstraint

startControlPoints_var = controlPoints_cvxpyVar[:,:d]
endControlPoints_var = controlPoints_cvxpyVar[:,-(d):]


#creates the equality contraints
startEqualityConstraint = [startControlPoints_var == initialControlPoints]
endEqualityConstraint = [endControlPoints_var == finalControlPoints]

controlPoints_constraints += startEqualityConstraint
controlPoints_constraints += endEqualityConstraint




#creates the objective function
velocityControlPoints_cp = controlPoints_cvxpyVar[:,0:-1] - controlPoints_cvxpyVar[:,1:]

#creates the objective function to adjust the control points for this thing.
minimizeLength_objectiveFunction = cp.Minimize(cp.sum(cp.norm(velocityControlPoints_cp, axis=1)))

prob = cp.Problem(minimizeLength_objectiveFunction, controlPoints_constraints)


#calls the function to solve the problem 
prob.solve(solver=cp.CLARABEL)

print("Status: ", prob.status)
print("objective: ", prob.value)
print("Optimized central control points: ", controlPoints_cvxpyVar.value)


controlPointsSolutions = controlPoints_cvxpyVar.value



bspline = BsplineEvaluation(control_points=controlPointsSolutions,
                            order = d,
                            start_time=0.0)

bsplinePoints, timePoints = bspline.get_spline_data(num_data_points_per_interval=100)

x_bsplinePoints = bsplinePoints[0,:]
y_bsplinePoints = bsplinePoints[1,:]


controlPoints_x = controlPointsSolutions[0,:]
controlPoints_y = controlPointsSolutions[1,:]


plt.figure(0)
plt.plot(x_bsplinePoints, y_bsplinePoints)
plt.scatter(controlPoints_x, controlPoints_y)
plt.xlim([0.0, 3.0])
plt.ylim([0.0, 6.0])
plt.show()


potato = 0