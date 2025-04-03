

import os, sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from bsplinegenerator.bsplines import BsplineEvaluation
from IPython.display import display

sys.path.insert(0,os.fspath(Path(__file__).parents[1]))
tempPath = sys.path


from eVTOL_BSplines.path_generation_helpers.conditions_helpers import conditions, conditionsList


#creates the conditions and the conditions list to test to see if they work


condition1 = conditions(numDerivatives=4, dimension=2)

#sets the conditions
condition1.setPosition(pos=np.array([[1],[1]]))
condition1.setVelocity(vel=np.array([[2],[2]]))
condition1.setAccel(accel=np.array([[3],[3]]))
condition1.setJerk(jerk=np.array([[4],[4]]))
condition1.setSnap(snap=np.array([[5],[5]]))


matrix = condition1.getConditionsMatrix()


condition2 = conditions(numDerivatives=4, dimension=2)

condition2.setPosition(pos=np.array([[11],[11]]))
condition2.setVelocity(vel=np.array([[22],[22]]))
condition2.setAccel(accel=np.array([[33],[33]]))
condition2.setJerk(jerk=np.array([[44],[44]]))
condition2.setSnap(snap=np.array([[55],[55]]))



#adds the condition to the matrix
list = conditionsList()
list.addCondition(condition1)
list.addCondition(condition2)



allConditions = list.getAllConditions()


tempCondition1 = list.getCondition(index=0)
tempCondition2 = list.getCondition(index=1)




potato = 0