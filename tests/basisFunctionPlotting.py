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

numSamplesPerSection = 100

sampler = basisFunctionSampler(degree=1,numSegmentsPerSection=numSamplesPerSection)

functionSampledList_init, initTimeList = sampler.getSampledData()

timeShift = 1.0

timeListsSections = []

#creates the knot vector
knotVector = [-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9]

pointsTime = knotVector[1:]

controlPoints1 = [1,2,-1,-0.5,0.5,1.5,3.0,5.0,4.5,-3.0]

controlPoints2 = deepcopy(controlPoints1)
#modifies one control point to hilight change
controlPoints2[5] = -2.0

scaledList1 = []
scaledList2 = []

timeListComplete = []


fig, ax = plt.subplots(2,3, sharex='col',figsize=(10, 6))

#iterates over all but the last item
for knot, ctrl1, ctrl2 in zip(knotVector[:-1], controlPoints1, controlPoints2):

    tempTimeList = [tempTime + knot for tempTime in initTimeList]

    #creates the scaled function for both sets of control points
    scaledFunction_1 = [ctrl1*tempSample for tempSample in functionSampledList_init]
    scaledFunction_2 = [ctrl2*tempSample for tempSample in functionSampledList_init]

    scaledList1.append(scaledFunction_1)
    scaledList2.append(scaledFunction_2)

    timeListsSections.append(tempTimeList)
    ax[0,0].plot(tempTimeList, functionSampledList_init, linewidth=4)
    ax[0,1].plot(tempTimeList, scaledFunction_1, linewidth=4)
    ax[0,2].plot(tempTimeList, scaledFunction_2, linewidth=4)

numSections = len(timeListsSections)

completeTimeList = []

for i, tempTimeList in enumerate(timeListsSections):
    
    #case we are a the beginning
    if i == 0:
        completeTimeList.extend(tempTimeList)
        testPoint = 0
    else:
        tempSection = tempTimeList[numSamplesPerSection:]
        completeTimeList.extend(tempSection)
        testPoint = 0


numTotalSamples = len(completeTimeList)

unmodifiedFunction = [0.0]*numTotalSamples
ctrl1_Function = [0.0]*numTotalSamples
ctrl2_Function = [0.0]*numTotalSamples

for i, (ctrl1Section, ctrl2Section) in enumerate(zip(scaledList1, scaledList2)):

    startIndex = i*numSamplesPerSection

    for j, (unmodVal, ctrl1Val, ctrl2Val) in enumerate(zip(functionSampledList_init, ctrl1Section, ctrl2Section)):

        currentIndex = startIndex + j
        unmodifiedFunction[currentIndex] = unmodifiedFunction[currentIndex] + unmodVal
        ctrl1_Function[currentIndex] = ctrl1_Function[currentIndex] + ctrl1Val
        ctrl2_Function[currentIndex] = ctrl2_Function[currentIndex] + ctrl2Val


testPoint = 0


#scatter plots the control points
ax[0,1].scatter(pointsTime, controlPoints1, label='Control Points 1', color='r')
ax[0,2].scatter(pointsTime, controlPoints2, label='Control Points 2', color='r')

ax[0,0].grid(True)
ax[0,0].legend()
ax[0,0].tick_params(axis='both', labelsize=14)
ax[0,0].set_title("Shifted Unscaled Basis Functions Degree 1")


ax[0,1].grid(True)
ax[0,1].legend()
ax[0,1].tick_params(axis='both', labelsize=14)
ax[0,1].set_title("Shifted Scaled Basis Functions Set 1")


ax[0,2].grid(True)
ax[0,2].legend()
ax[0,2].tick_params(axis='both', labelsize=14)
ax[0,2].set_title("Shifted Scaled Basis Functions Set 2")

#plots the summation plots
ax[1,0].plot(completeTimeList, unmodifiedFunction, color='green', linewidth=4)
ax[1,0].legend()
ax[1,0].grid(True)
ax[1,0].tick_params(axis='both', labelsize=14)
ax[1,0].set_title("Shifted Unscaled Basis Sum Degree 1")

ax[1,1].plot(completeTimeList, ctrl1_Function, color='green', linewidth=4)
ax[1,1].scatter(pointsTime, controlPoints1, label='Control Points 1', color='r')
ax[1,1].legend()
ax[1,1].grid(True)
ax[1,1].tick_params(axis='both', labelsize=14)
ax[1,1].set_title("Control Points 1 Sum")


ax[1,2].plot(completeTimeList, ctrl2_Function, color='green', linewidth=4)
ax[1,2].scatter(pointsTime, controlPoints2, label='Control Points 2', color='r')
ax[1,2].legend()
ax[1,2].grid(True)
ax[1,2].tick_params(axis='both', labelsize=14)
ax[1,2].set_title("Control Points 2 Sum")

plt.suptitle("Degree 1 Basis Functions with Examples")

plt.show()


testPoint = 0


