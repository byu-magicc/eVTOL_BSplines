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



d = 5
M = 10


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


V = V_t.T

Sigma_inv = np.linalg.inv(np.diag(Sigma))

U1_t = U1.T

print("V")
display(sp.Matrix(np.round(V,decimals=3)))
print("Sigma inv")
display(sp.Matrix(np.round(Sigma_inv,decimals=3)))
print("U_t")
display(sp.Matrix(np.round(U1_t,decimals=3)))


#gets the inverse portion
inverse_portion = V @ Sigma_inv @ U1_t
print("Inverse Portion")
display(sp.Matrix(np.round(inverse_portion, decimals=3)))

#print("Inverse Portion")
#display(inverse_portion)



U2_mult = U2 @ np.transpose(U2)
print("U_2 U_2^t")
display(sp.Matrix(np.round(U2_mult, decimals=3)))

potato = 0





# %%

#gets B matrix, both normal and inverse
B_forward = B

print("B forward")
display(sp.Matrix(np.round(B_forward, decimals=3)))

#gets the B reverse portion
B_reverse = inverse_portion
print("B reverse")
display(sp.Matrix(np.round(B_reverse, decimals=3)))



#gets the B hat for the same parameters
B_hat, B_hat_inv = B_hat_B_hat_inv(degree=d, alpha=1.0, M=M)

print("B hat")
display(sp.Matrix(np.round(B_hat,decimals=3)))
print("B hat inv")
display(sp.Matrix(np.round(B_hat_inv,decimals=3)))





#creates the function to create the w matrix
w_matrix_creator = create_W_Matrix(d=d, integratorFileName = "lookUpTables/degree_5_integrations.npz")

rho = np.ones((d))

#gets the full W matrix
W_matrix = w_matrix_creator.W_d_M_rho(d=d, M=M, rho=rho)

print("W_matrix")
display(sp.Matrix(np.round(W_matrix, decimals=3)))
print(sp.latex(sp.Matrix(np.round(W_matrix, decimals=3))))
#gets the W matrix partition
W_matrix_partitioned = w_matrix_creator.W_partitioned(d=d, M=M, rho=rho)

for i in range(3):
    for j in range(3):
        temp_W = W_matrix_partitioned[i][j]
        display(sp.Matrix(np.round(temp_W, decimals=3)))

#gets the W matrix inverse
W_matrix_inv = np.linalg.inv(W_matrix)
print("W inv")
display(sp.Matrix(np.round(W_matrix_inv, decimals=3)))
print(sp.latex(sp.Matrix(np.round(W_matrix_inv, decimals=3))))

#gets the center section of the inverse
W_inv_center = W_matrix_inv[d:(M),d:(M)]
display(sp.Matrix(np.round(W_inv_center, decimals=3)))
print(sp.latex(sp.Matrix(np.round(W_inv_center, decimals=3))))

#gets the inverse of the center of the original W matrix
W_center_inv = np.linalg.inv(W_matrix_partitioned[1][1])

display(sp.Matrix(np.round(W_center_inv, decimals=3)))
print(sp.latex(sp.Matrix(np.round(W_center_inv, decimals=3))))
# %%

#this section gets the U2 Transpose U2
U2_t_U2 = np.transpose(U2) @ U2

print("U2 transpose U2")
display(sp.Matrix(np.round(U2_t_U2, decimals=3)))




# %%
