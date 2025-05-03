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



k = 3
M = 10


#gets the S matrix through a test

S_k_m = W_matrix_Creation.S_k_M(k=k,
                                M=M)


S_k_m_actual = integrate_b_bT(degree=k,
                              l=0,
                              M=M)

error = S_k_m_actual - S_k_m


#gets the sum of the error

display(sp.Matrix(np.round(S_k_m, decimals=3)))

display(sp.Matrix(np.round(S_k_m_actual, decimals=3)))

display(sp.Matrix(error))

A_temp = (W_matrix_Creation.integration_Matrix_List)[3]

diagonal = W_matrix_Creation.getIMatrixDiagonal(I=A_temp, index=0, isColumn=False)

section = W_matrix_Creation.getDiagonalSection(removingLeft=False, numRemovals=0, diagonal=diagonal)

W_matrix_Creation.sumDiagonalTopRow(A=A_temp, index=0)

W_matrix_Creation.getDiagonalLength(A=A_temp, index = 0)


hurricane = 0


# %%
