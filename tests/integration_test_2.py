#%%
import os, sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from bsplinegenerator.bsplines import BsplineEvaluation
from IPython.display import display

from IPython.display import display

import scipy.integrate as integrate

sys.path.insert(0,os.fspath(Path(__file__).parents[1]))
tempPath = sys.path


from eVTOL_BSplines.path_generation_helpers.conditions_helpers import conditions, conditionsList

from eVTOL_BSplines.path_generation_helpers.matrix_helpers import S_M_d_l, get_W_d_l_M, integrate_b_bT


#sets the degree, l and M
degree = 3
l = 2
M = 5

#gets the S Matrix

S = S_M_d_l(M=M, degree=degree, l=(l-1))

#gets the integration matrix
integratedMatrix = integrate_b_bT(degree=degree, l=l, M=M)

display(S)
display(np.shape(S))

display(integratedMatrix)
display(np.shape(integratedMatrix))
# %%
