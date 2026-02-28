import os, sys
from pathlib import Path
import numpy as np

import matplotlib.pyplot as plt
from sympy.printing.pretty.stringpict import line_width
from copy import deepcopy

sys.path.insert(0,os.fspath(Path(__file__).parents[1]))
tempPath = sys.path

from bsplinegenerator.bsplines import BsplineEvaluation
from eVTOL_BSplines.path_generation_helpers.basis_function_helpers import basisFunctionSampler


sampler = basisFunctionSampler(degree=1,numSegmentsPerSection=100)

functionSampledList_init, initTimeList = sampler.getSampledData()

timeShift = 1.0

timeListsComplete = []

#creates the knot vector
knotVector = [-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9]

pointsTime = knotVector[1:]

controlPoints1 = [1,2,-1,-0.5,0.5,1.5,3.0,5.0,4.5,-3.0]

controlPoints2 = deepcopy(controlPoints1)
#modifies one control point to hilight change
controlPoints2[5] = -2.0

scaledList1 = []
scaledList2 = []

#unmodified basis functions
fig, ax1 = plt.subplots(1, 1, sharex=True)
#figure for first set of control points
fig, ax2 = plt.subplots(1, 1, sharex=True)
#figure for second set of control points
fig, ax3 = plt.subplots(1, 1, sharex=True)

#iterates over all but the last item
for knot, ctrl1, ctrl2 in zip(knotVector[:-1], controlPoints1, controlPoints2):

    tempTimeList = [tempTime + knot for tempTime in initTimeList]

    #creates the scaled function for both sets of control points
    scaledFunction_1 = [ctrl1*tempSample for tempSample in functionSampledList_init]
    scaledFunction_2 = [ctrl2*tempSample for tempSample in functionSampledList_init]

    scaledList1.append(scaledFunction_1)
    scaledList2.append(scaledFunction_2)

    timeListsComplete.append(tempTimeList)
    ax1.plot(tempTimeList, functionSampledList_init, linewidth=4)
    ax2.plot(tempTimeList, scaledFunction_1, linewidth=4)
    ax3.plot(tempTimeList, scaledFunction_2, linewidth=4)

    #scatter plots the control points
    ax2.scatter(pointsTime, controlPoints1, label='Control Points 1')
    ax3.scatter(pointsTime, controlPoints2, label='Control Points 2')


ax1.grid(True)
ax1.legend()
ax1.set_ylim((-0.5, 1.5))
ax1.tick_params(axis='both', labelsize=14)
ax1.set_title("Shifted Unscaled Basis Functions Degree 1")


ax2.grid(True)
ax2.legend()
ax2.tick_params(axis='both', labelsize=14)
ax2.set_title("Shifted Scaled Basis Functions Set 1")


ax3.grid(True)
ax3.legend()
ax3.tick_params(axis='both', labelsize=14)
ax3.set_title("Shifted Scaled Basis Functions Set 2")


plt.show()


testPoint = 0


