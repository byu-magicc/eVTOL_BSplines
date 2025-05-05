#%%
import os, sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from bsplinegenerator.bsplines import BsplineEvaluation
from IPython.display import display
import sympy as sp

from IPython.display import display

import scipy.integrate as integrate

sys.path.insert(0,os.fspath(Path(__file__).parents[1]))
tempPath = sys.path


from eVTOL_BSplines.path_generation_helpers.conditions_helpers import conditions, conditionsList

from eVTOL_BSplines.path_generation_helpers.matrix_helpers import b_d_M_t_vector, integrate_b_bT

from eVTOL_BSplines.path_generation_helpers.matrix_generators_efficient import create_W_Matrix

degree = 5

W_matrix_Creation = create_W_Matrix(d=degree,
                                      integratorFileName = "lookUpTables/degree_5_integrations.npz")


Num_m = 50
Num_d = 5


k = 3
M = 9


#let's do a brute force method
def is_toeplitz(matrix: np.ndarray):
    matrix = np.array(matrix)
    return all(np.diag(matrix, k).size == 1 or np.all(np.diag(matrix, k) == np.diag(matrix, k)[0])
               for k in range(-matrix.shape[0] + 1, matrix.shape[1]))





#creates list to store whether each test matrix is toeplitz
toeplitzChecks = []

for i in range(10, Num_m):
    for j in range(1, Num_d+1):
        
        #iterating over all of these M and d values, we check whether or not the resulting matrix is toeplitz
        W_temp = W_matrix_Creation.W_d_l_M(d=j,l=1,M=i)

        #gets the partition of W_temp
        W_temp_partition = W_matrix_Creation.getWMatrixPartition(W=W_temp, d=j, M=i)
        W_temp_center = (W_temp_partition[1])[1]
        tempStatus = is_toeplitz(matrix=W_temp_center)

        toeplitzChecks.append(tempStatus)


#checks if the list contains at least one false
if not all(toeplitzChecks):
    print("one element is false")
else:
    print("all elements are true")
        
hurricane = 0




# %%
