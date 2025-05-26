import os, sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from bsplinegenerator.bsplines import BsplineEvaluation
from IPython.display import display

sys.path.insert(0,os.fspath(Path(__file__).parents[1]))
tempPath = sys.path


from eVTOL_BSplines.path_generation_helpers.lookUpHelpers import lookUpTableGenerator

from eVTOL_BSplines.path_generation_helpers.matrix_helpers import b_d_M_t_vector


b0 = b_d_M_t_vector(time=0.0, degree=1, alpha=1.0, M=4)
b4 = b_d_M_t_vector(time=4.0, degree=1, alpha=1.0, M=4)



#creates the main function to implement the look up table generator

generator = lookUpTableGenerator(M_maximum=100, highestDegree=5)

#calls the function to generate the lookup tables
generator.generateSLookupTables()
generator.generateWLookupTables()
generator.generateBEndTables()
generator.generateYZLookupTables()


taco = 0
