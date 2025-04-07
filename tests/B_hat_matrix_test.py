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



#imports the bspline
from eVTOL_BSplines.path_generation_helpers.matrix_helpers import B_hat_B_hat_inv_simplified, B_hat_B_hat_inv



#sets the variables
M = 5
d = 3

#gets the old fashioned B had and inverse using the old fashioned method
B_hat, B_hat_inv = B_hat_B_hat_inv(degree=d, alpha=1.0, M=M)



B_hat_simp, B_hat_inv_simp = B_hat_B_hat_inv_simplified(degree=d, alpha=1.0)




#gets the error
error = B_hat - B_hat_simp

error_inv = B_hat_inv - B_hat_inv_simp


#displays the B_hat
display("B hat")

display(B_hat_simp)

display("B hat inv")

display(B_hat_inv_simp)


display(error)
display(error_inv)
# %%
