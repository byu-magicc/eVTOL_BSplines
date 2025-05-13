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
from eVTOL_BSplines.path_generation_helpers.matrix_helpers import B_M_matrix, B_hat_B_hat_inv



d = 3
M = 9



#gets the B_M matrix concatenated together
B_0 = B_M_matrix(time=0, degree=d, alpha=1.0, M=M)
B_M = B_M_matrix(time=M, degree=d, alpha=1.0, M=M)

#puts them together
B = np.concatenate((B_0, B_M), axis=1)
print("B")
display(sp.Matrix(B))



#gets the SVD of the B matrix
U, Sigma, V_t = np.linalg.svd(B)



#gets the U1 and U2
U1 = U[:,:(2*d)]
U2 = U[:,(2*d):]
print("U_1")
display(sp.Matrix(np.round(U1,decimals=3)))
print("U_2")
display(sp.Matrix(np.round(U2,decimals=3)))

print("Sigma")
display(sp.Matrix(np.round(Sigma,decimals=3)))

print("V_t")
display(sp.Matrix(np.round(V_t,decimals=3)))

#creates the function to create the w matrix
w_matrix_creator = create_W_Matrix(d=d, integratorFileName = "lookUpTables/degree_5_integrations.npz")

rho = np.ones((d))

#obtains the W matrix
W_matrix = w_matrix_creator.W_d_M_rho(d=d, M=M, rho=rho)

print("W matrix")
display(sp.Matrix(np.round(W_matrix,decimals=3)))



#gets the inverted matrix
inversionInput = np.transpose(U2) @ W_matrix @ U2

print("input to inversion")
display(sp.Matrix(np.round(inversionInput,decimals=3)))

inversionOutput = np.linalg.inv(inversionInput)
print("output from inversion")
display(sp.Matrix(np.round(inversionOutput,decimals=3)))


# %%
