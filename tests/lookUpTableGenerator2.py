import os, sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from bsplinegenerator.bsplines import BsplineEvaluation
from IPython.display import display

sys.path.insert(0,os.fspath(Path(__file__).parents[1]))
tempPath = sys.path

from eVTOL_BSplines.path_generation_helpers.lookUpHelpers_2 import lookUpTablesGenerator, YZGeneratorReader

#sets the maximum M and degree
M_max = 100
d_max = 5


#creates the generator 
generator = lookUpTablesGenerator(M_maximum=M_max, 
                                  highestD=d_max)

yzGen = YZGeneratorReader()

#calls the function to generate the tables associated with the B end tables
generator.generateBEndTables()

#calls the function to generate the S matrix file
generator.generateSLookupTables()


#calls the function to generate the W matrix file
generator.generateWLookupTables()


#calls the generator function for the YZ lookup tables
yzGen.generateYZLookupTables()


taco = 0