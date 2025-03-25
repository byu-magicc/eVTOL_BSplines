#%%

import os, sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from bsplinegenerator.bsplines import BsplineEvaluation
from IPython.display import display

sys.path.insert(0,os.fspath(Path(__file__).parents[1]))
tempPath = sys.path

from eVTOL_BSplines.path_generation_helpers.matrix_helpers import uniform_cox_de_Boor_basis_function, uniform_knot_point_generator, uniform_cox_de_boor_basis_function_table
from eVTOL_BSplines.path_generation_helpers.matrix_helpers import uniform_basis_function_evaluation
from eVTOL_BSplines.path_generation_helpers.matrix_helpers import b_d_M_t_vector, B_d_M_t_matrix

#gets the b vector
b_d_M = b_d_M_t_vector(time=0.0,
                       degree=2,
                       M=6)

display(b_d_M)


#gets the B matrix
B_d_M, B_hat_d_M = B_d_M_t_matrix(time=1.0,
                       degree=3,
                       M=5)

display(B_d_M)
display(B_hat_d_M)

# %%
