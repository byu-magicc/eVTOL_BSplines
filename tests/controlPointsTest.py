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


from eVTOL_BSplines.path_generation_helpers.matrix_generators_efficient import create_W_Matrix

degree = 5
M=50

rho = np.array([[1.0,1.0,1.0,1.0,1.0]])

W_matrix_Creation = create_W_Matrix(d=degree,
                                      integratorFileName = "lookUpTables/degree_5_integrations.npz")


#gets the partitioned matrix
W_partitioned = W_matrix_Creation.W_partitioned(d=degree,M=50,rho=rho)



potato = 0