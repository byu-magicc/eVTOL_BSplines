#does the same thing as the control points matrix file, but accomplishes this
#using lookup tables instead of recursions


import numpy as np
import os, sys
from pathlib import Path
#imports the things I have 
from eVTOL_BSplines.path_generation_helpers.matrix_helpers import D_d_l_M, B_init_final, B_init_final_svd

