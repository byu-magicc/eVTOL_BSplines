import os, sys
from pathlib import Path
import numpy as np

sys.path.insert(0,os.fspath(Path(__file__).parents[1]))
tempPath = sys.path


from eVTOL_BSplines.path_generation_helpers.staticFlightPath import staticFlightPath
from eVTOL_BSplines.path_generation_helpers.localControlPoints import getLocalControlPoints




#creates the initial conditionss

startPosition = np.array([[0.0],
                          [0.0],
                          [0.0]])

startVelocity = np.array([[25.0],
                          [0.0],
                          [0.0]])

startAcceleration = np.array([[0.0],
                              [0.0],
                              [0.0]])



endPosition = np.array([[100.0],
                        [50.0],
                        [0.0]])

endVelocity = np.array([[0.0],
                        [25.0],
                        [0.0]])

endAccleration = startAcceleration


startConditions = [startPosition, startVelocity, startAcceleration]
endConditions = [endPosition, endVelocity, endAccleration]

#creates the rho list
rho = np.array([[1.0],
                [1.0],
                [1.0]])


staticGen = staticFlightPath()

#gets the control points
controlPoints = staticGen.getControlPoints(initialConditions=startConditions,
                                           finalConditions=endConditions,
                                           rho=rho,
                                           numDimensions=3,
                                           d=3,
                                           M=10)



potato = 0