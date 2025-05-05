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



#gets the W matrix
W = W_matrix_Creation.W_d_l_M(d=k, l=1, M=M)
display(sp.Matrix(np.round(W,decimals=3)))

#gets the S matrix through a test

S_k_m, numSectionsMatrix = W_matrix_Creation.S_k_M(k=k,
                                M=M)


S_k_m_actual = integrate_b_bT(degree=k,
                              l=0,
                              M=M)

error = S_k_m_actual - S_k_m

W_partition = W_matrix_Creation.getWMatrixPartition(W=W, d=k, M=M)

temp = W_partition[1]

display(sp.Matrix(np.round((W_partition[1])[1], decimals=3)))


#let's do a brute force method


def is_toeplitz(matrix):
    matrix = np.array(matrix)
    return all(np.diag(matrix, k).size == 1 or np.all(np.diag(matrix, k) == np.diag(matrix, k)[0])
               for k in range(-matrix.shape[0] + 1, matrix.shape[1]))




hurricane = 0




# %%
