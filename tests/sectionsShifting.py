#%%
import sympy as sp
from IPython.display import display

t_1 = sp.symbols('t_1')
t_2 = sp.symbols('t_2')
t_3 = sp.symbols('t_3')



t = sp.symbols('t')
#degree 3
section_1 = sp.expand((-1/2)*(t+1)**3 + 2*(t+1)**2 - 2*(t+1) + (2/3))
display(section_1)

section_2 = sp.expand((1/2)*(t+2)**3 - 4*(t+2)**2 + 10*(t+2) - (22/3))
display(section_2)

# %%
