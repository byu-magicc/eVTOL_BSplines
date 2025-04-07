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
import sympy as sp

sys.path.insert(0,os.fspath(Path(__file__).parents[1]))
tempPath = sys.path


from eVTOL_BSplines.path_generation_helpers.conditions_helpers import conditions, conditionsList

from eVTOL_BSplines.path_generation_helpers.matrix_helpers import get_W_d_l_M, get_W_d_M_rho, get_W_partitioned


#tests the W matrices
M = 5
degree = 2
l = degree - 1


#creates the rho matrix
rho = np.array([[0],[1]])

individual_w = get_W_d_l_M(degree=degree, l=l, M=M)
total_W = get_W_d_M_rho(degree=degree, L=l, M=M, rho=rho)
partitioned_W = get_W_partitioned(degree=degree, M=M, L=l, rho=rho)




total_W = np.round(total_W, decimals=3)

#iterates over every submatrix in the matrix list from the partitioned W
for row in partitioned_W:
    for i, arr in enumerate(row):
        row[i] = np.round(arr, 3)


total_W_sym = sp.Matrix(total_W)
partitioned_W_sym = sp.Matrix(partitioned_W)

display(total_W_sym)
display(partitioned_W_sym)


carrot = 0
# %%
