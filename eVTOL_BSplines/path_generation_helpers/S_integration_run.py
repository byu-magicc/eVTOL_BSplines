#%%

from IPython.display import display
import sympy as sp
import numpy as np
from S_integration import S_integration_naieve




integrator = S_integration_naieve(degree=5, M=7)




roundedMatrix = np.round(integrator.S_m_d, decimals=5)
display(sp.Matrix(roundedMatrix))


potat = 0

# %%
