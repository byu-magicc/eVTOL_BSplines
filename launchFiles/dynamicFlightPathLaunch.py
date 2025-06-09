#implements the launch file for the dynamic flight path creator
import os, sys
from pathlib import Path
import numpy as np
import time

sys.path.insert(0,os.fspath(Path(__file__).parents[1]))
tempPath = sys.path


numDimensions = 2
numConditions = 3

degree = 3

#creates the rho array
rho = np.array([1.0, 1.0, 1.0])

M=15

from eVTOL_BSplines.path_generation_helpers.dynamicallyAdjustingFlightPath import dynamicFlightPath
from eVTOL_BSplines.path_generation_helpers.conditions import conditions
#imports the BSplineEvaluation function for this thing.
from bsplinegenerator.bsplines import BsplineEvaluation




flightConditions = conditions(dimension=numDimensions,
                              numConditions=numConditions)


#creates the conditions list
pos_init = np.array([[0.0], [0.0]])
vel_init = np.array([[1.0], [0.0]])
accel_init = np.array([[2.0], [0.0]])


conditionsList_init = [pos_init, vel_init, accel_init]


#creates the final conditions list
pos_final = np.array([[50.0],[10.0]])
vel_final = np.array([[0.0],[5.0]])
accel_final = np.array([[0.0],[1.0]])


conditionsList_final = [pos_final, vel_final, accel_final]


flightGen = dynamicFlightPath(initialConditionsMain=conditionsList_init,
                                finalConditionsMain=conditionsList_final,
                                numDimensions=numDimensions,
                                d=degree,
                                M=M)

#and let's time how long each iteration takes


#calls the function to get the inverted portion
startTime = time.time()
CtrlPnts = flightGen.getControlPoints(rho=rho)
endTime = time.time()
timeDifference = endTime - startTime
print(timeDifference)

#creates the bspline class, which is useful to use here
outputBSpline = BsplineEvaluation(control_points=CtrlPnts,
                                  order=degree,
                                  start_time=0.0)

outputBSpline.plot_spline(num_data_points_per_interval=100)

#gets the minvo control points
MinvoCtrlPts = outputBSpline.get_minvo_control_points()

outputBSpline.plot_minvo_curves(num_data_points_per_interval=100)

potato = 0