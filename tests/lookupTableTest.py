#implements the test file for the lookup tables

import os, sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from bsplinegenerator.bsplines import BsplineEvaluation
from IPython.display import display

sys.path.insert(0,os.fspath(Path(__file__).parents[1]))
tempPath = sys.path

from eVTOL_BSplines.path_generation_helpers.control_points_lookup import controlPointsGeneratorLookup


M = 5

d = 2


dimension = 2


#creates the generator
generator = controlPointsGeneratorLookup(M=M,
                                         degree=d,
                                         dimension=dimension)




#calls the Y summation function

Y_sum = generator.sumYArrays(rho=np.array([1,1]))




Z_sum = generator.sumZArrays(rho=np.array([1,1]))




ashes = 0