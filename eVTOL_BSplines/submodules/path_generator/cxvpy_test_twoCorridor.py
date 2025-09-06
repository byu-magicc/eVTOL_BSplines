#this file implements a two corridor thing for the optimal control point

#this file is my first preliminary test to see how cvxpy actually works
import cvxpy as cp
import numpy as np
from numpy import cos, sin


from bsplinegenerator.bsplines import BsplineEvaluation

import matplotlib.pyplot as plt

M = 10

d = 3



#creates the bounds for the first convex hull
vertices_1 = np.array([[1.0, 1.0, 2.0, 2.0],
                       [1.0, 6.0, 6.0, 1.0]])


#creates the A1 vector
normals_1 = np.array([[-1.0, 0.0, 1.0, 0.0],
                      [0.0, 1.0, 0.0, -1.0]])

A1 = normals_1.T


b1 = np.array([[-1.0],
               [6.0],
               [2.0],
               [-1.0]])


#creates the stuff for the second corridor
vertices_2 = np.array([[1.0, 6.0, 6.0, 1.0],
                       [6.0, 6.0, 5.0, 5.0]])



normals_2 = np.array([[0.0, 1.0, 0.0, -1.0],
                      [1.0, 0.0, -1.0, 0.0]])



A2 = normals_2.T

b2 = np.array([[6.0],
               [6.0],
               [-5.0],
               [-1.0]])



initialControlPoints = np.array([[1.5, 1.5, 1.5],
                                 [1.25, 1.5, 1.75]])

finalControlPoints = np.array([[5.25, 5.5, 5.75],
                               [5.5, 5.5, 5.5]])



#gets the total number of control points
numControlPoints = 2*M + d


#creates the control points variables
controlPoints_cvxpyVar = cp.Variable((2, numControlPoints))


#creates the two partitions of the control points
controlPoints_cvxpyVar_1 = controlPoints_cvxpyVar[:,:(M+d)]

controlPoints_cvxpyVar_2 = controlPoints_cvxpyVar[:,M:]


controlPoints_constraints = []


boundaryConstraint_1 = [A1 @ controlPoints_cvxpyVar_1 <= b1]
boundaryConstraint_2 = [A2 @ controlPoints_cvxpyVar_2 <= b2]

controlPoints_constraints += boundaryConstraint_1
controlPoints_constraints += boundaryConstraint_2


startControlPoints_var = controlPoints_cvxpyVar[:,:d]
endControlPoints_var = controlPoints_cvxpyVar[:,(-d):]


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
plt.scatter(controlPoints_x, controlPoints_y, color='orange')
plt.xlim([0.0, 7.0])
plt.ylim([0.0, 7.0])
plt.show()


potato = 0