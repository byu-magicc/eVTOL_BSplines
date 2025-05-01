#%%


import numpy as np
import sympy as sp

from matrix_helpers import *

from matrix_generators_efficient import *
from IPython.display import display

M = 6
d = 4


l = 2

D_result = D_d_l_M(d=d, l=l, M=M)




potato = 0



a0, a1, a2, a3, a4 = sp.symbols('a_0, a_1, a_2, a_3, a_4')
b0, b1, b2, b3, b4 = sp.symbols('b_0, b_1, b_2, b_3, b_4')


A = sp.Matrix([[a0, a1, a2, a3, a4],
               [a1, a0, a1, a2, a3],
               [a2, a1, a0, a1, a2],
               [a3, a2, a1, a0, a1],
               [a4, a3, a2, a1, a0]])

B = sp.Matrix([[b0, b1, b2, b3, b4],
               [b1, b0, b1, b2, b3],
               [b2, b1, b0, b1, b2],
               [b3, b2, b1, b0, b1],
               [b4, b3, b2, b1, b0]])



C = A*B
display(C)
display(B*A)


# %%
