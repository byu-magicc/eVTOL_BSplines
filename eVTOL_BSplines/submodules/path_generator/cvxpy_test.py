#this file is my first preliminary test to see how cvxpy actually works
import cvxpy as cp
import numpy as np
from numpy import cos, sin


#sets the number of points
numPoints = 4
#sets the number of dimensions
numDimensions = 2
#creates the cp variable from that number of points
points = cp.Variable((numPoints, numDimensions))

#gets the center of the rotated rectangle
rectangleCenter = np.array([10.0, 10.0])
#sets the width
rectangleWidth = 3.0
#sets the height of the rectangle
rectangleHeight = 1.5

#sets the rotation angle
theta = np.pi/6

#gets the rotation matrix from the corridor frame to the world frame
Rotation_corridorToWorld = np.array([[cos(theta), -sin(theta)],
                                     [sin(theta),  cos(theta)]])

#gets the rotation matrix from the world frame to the corridor frame
Rotation_worldToCorridor = Rotation_corridorToWorld.T


#goes through and gets the constraings
constraints = []
for i in range(numPoints):
    #gets the current point variable (the independent variable set by the cp class)
    currentVariablePoint = points[i]
    #gets the vector from the center to the current point variable in the world frame
    centerToCurrentVariablePoint_world = currentVariablePoint - rectangleCenter
    #gets the center to current variable point in the corridor frame
    centerToCurrentVariablePoint_corridor = Rotation_worldToCorridor @ centerToCurrentVariablePoint_world

    #gets the current constraint for x and y (in the corridor frame. )
    currentConstraint = [centerToCurrentVariablePoint_corridor[0] <= rectangleWidth/2.0,
                         centerToCurrentVariablePoint_corridor[0] >= -rectangleWidth/2.0,
                         centerToCurrentVariablePoint_corridor[1] <= rectangleHeight/2.0,
                         centerToCurrentVariablePoint_corridor[1] >= -rectangleHeight/2.0]
    #adds the current constraint to the constraints list
    constraints += currentConstraint

#creates the objective function
objective = cp.Minimize(cp.sum(cp.norm(points[1:] - points[:-1], axis=1)))

#gets the problem
prob = cp.Problem(objective=objective,
                  constraints=constraints)

#gets the solve problem
prob.solve(solver=cp.ECOS)

print("Optimal points:\n", points.value)


potato = 0