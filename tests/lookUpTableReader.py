import os, sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from bsplinegenerator.bsplines import BsplineEvaluation
from IPython.display import display

sys.path.insert(0,os.fspath(Path(__file__).parents[1]))
tempPath = sys.path


from eVTOL_BSplines.path_generation_helpers.lookUpHelpers import lookUpTableReader


#creates the reader
reader = lookUpTableReader()
#reads from the S lookup table
S_table = reader.readSLookupTables()

#reads from the B lookup table
B_table = reader.readBLookupTable()


#gets the part of the B table I want to analyze
B_temp = B_table['d2_M20_B']
U1_temp = B_table['d2_M20_U1']
U2_temp = B_table['d2_M20_U2']
Sigma = B_table['d2_M20_Sigma']
Vt = B_table['d2_M20_Vt']



potato = 0