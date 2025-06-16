#This file implements the dynamic flight path launch file, which will be an animation 
#that updates regularly. In our case, I think it will just draw conditions from the
#original or initial flight path and just be shorter segments of that.

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


M = 15


from eVTOL_BSplines.path_generation_helpers.dynamicallyAdjustingFlightPath import dynamicFlightPath
from eVTOL_BSplines.path_generation_helpers.conditions import conditions
from bsplinegenerator.bsplines import BsplineEvaluation



#defines a simple sectioning function for the knot points
def sectionKnots(knots: np.ndarray,
                 M: int,
                 d: int):
    
    #gets the preSection
    preSection = knots[:d]

    #gets the main section
    mainSection = knots[d:(M+d+1)]

    #gets the post section
    postSection = knots[(M+d+1):]


    return preSection, mainSection, postSection



flightConditions = conditions(dimension=numDimensions,
                              numConditions=numConditions)


#creates the conditions list
pos_init = np.array([[0.0], [0.0]])
vel_init = np.array([[0.0], [5.0]])
accel_init = np.array([[0.0], [1.0]])


conditionsList_init = [pos_init, vel_init, accel_init]


#creates the final conditions list
pos_final = np.array([[50.0],[10.0]])
vel_final = np.array([[10.0],[0.0]])
accel_final = np.array([[0.1],[0.0]])


conditionsList_final = [pos_final, vel_final, accel_final]



#instantiates the flight generator
flightGen = dynamicFlightPath(initialConditionsMain=pos_init,
                              finalConditionsMain=conditionsList_final)


#gets the Control points for the initial and final state specified
CtrlPnts = flightGen.getCurrentControlPoints(current_M=M,
                                             rho=rho,
                                             currentInitialConditions=conditionsList_init,
                                             currentFinalConditions=conditionsList_final)


#creates the 
outputBSpline = BsplineEvaluation(control_points=CtrlPnts,
                                  order=degree,
                                  start_time=0.0)

#gets the knot points
knotPnts = outputBSpline.get_knot_points()

preKnots, mainKnots, postKnots = sectionKnots(knots=knotPnts,
                                              M=M,
                                              d=degree)

#creates the list of the full conditions list
fullConditionsList = []

#iterates over the main knots
for knot in mainKnots:

    #gets the spline at the knot point time
    tempPosition = outputBSpline.get_spline_at_time_t(time=knot)

    #gets the velocity
    tempVelocity = outputBSpline.get_derivative_at_time_t(time=knot,
                                                      derivative_order=1)
    
    #gets the acceleration
    tempAcceleration = outputBSpline.get_derivative_at_time_t(time=knot,
                                                          derivative_order=2)

    #puts them together into the temp conditions list
    tempConditionsList = [tempPosition, tempVelocity, tempAcceleration]

    #appends to the full conditions list
    fullConditionsList.append(tempConditionsList)

    tomato = 0


#creates the full control points list
fullControlPointsList = []

#creates the bspline list
fullBSplineList = []

#iterates over all the items in the conditions list, and gets the unique control point set for all of them
for i, tempCondition in enumerate(fullConditionsList):

    #gets the current M
    current_M = M - i

    #checks if we have gone beyond the number acceptable for control points
    #at which point we break from the animation
    if current_M <= (degree):
        break

    #using the conditions and the final condition, we get the appropriate update for the bspline
    tempControlPoints = flightGen.getCurrentControlPoints(current_M=current_M,
                                                          rho=rho,
                                                          currentInitialConditions=tempCondition,
                                                          currentFinalConditions=conditionsList_final)
    
    #saves the control points here
    fullControlPointsList.append(tempControlPoints)


    #generates the bspline
    tempBspline = BsplineEvaluation(control_points=tempControlPoints,
                                    order=degree,
                                    start_time=0.0)
    
    fullBSplineList.append(tempBspline)

    #plots the bspline
    tempBspline.plot_spline(num_data_points_per_interval=100)


#iterates over all of the control points as they get progressively smaller and smaller
#iterates over the prospective number of M's

#and gets the conditions from the spline at each of the control points
for temp_M in range(M, (degree - 1), -1):
    print(temp_M)