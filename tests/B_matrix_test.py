#I'm just using this file to test the B matrices generation


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
from eVTOL_BSplines.path_generation_helpers.matrix_helpers import get_W_partitioned


M = 15
d = 10
L = d - 1




potato = 0