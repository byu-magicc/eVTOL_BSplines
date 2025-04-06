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

from eVTOL_BSplines.path_generation_helpers.matrix_helpers import get_W_d_l_M, get_W_d_M_rho


#tests the W matrices
M = 5
degree = 2
l = 1


#creates the rho matrix
rho = np.array([[0],[1]])

individual_w = get_W_d_l_M(degree=degree, l=l, M=M)
total_W = get_W_d_M_rho(degree=degree, L=l, M=M, rho=rho)

display(individual_w)
display(total_W)
display(individual_w - total_W)

carrot = 0
# %%
