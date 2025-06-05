import os, sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from bsplinegenerator.bsplines import BsplineEvaluation
from IPython.display import display

sys.path.insert(0,os.fspath(Path(__file__).parents[1]))
tempPath = sys.path


from eVTOL_BSplines.path_generation_helpers.lookUpHelpers_2 import lookUpTableReader


#instantiates the reader
reader = lookUpTableReader()


B, U1, U2, Sigma, Vt, Pseudoinverse = reader.loadBLookupTables()

B_temp = reader.getIndividualB(M=10, d=2)
U1_temp = reader.getIndividualU1(M=10, d=2)
U2_temp = reader.getIndividualU2(M=10, d=2)
Sigma_temp = reader.getIndividualSigma(M=10, d=2)
Vt_temp = reader.getIndividualVt(M=10, d=2)
Pseudoinverse_temp = reader.getIndividualPseudoinverse(M=10, d=2)



tomato = 0