import os, sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from bsplinegenerator.bsplines import BsplineEvaluation
from IPython.display import display

sys.path.insert(0,os.fspath(Path(__file__).parents[1]))
tempPath = sys.path


from eVTOL_BSplines.path_generation_helpers.lookUpHelpers_2 import lookUpTableReader, YZGeneratorReader


#instantiates the reader
reader = lookUpTableReader()

YZReader = YZGeneratorReader()

'''
B, U1, U2, Sigma, Vt, Pseudoinverse = reader.loadBLookupTables()

B_temp = reader.getIndividualB(M=10, d=2)
U1_temp = reader.getIndividualU1(M=10, d=2)
U2_temp = reader.getIndividualU2(M=10, d=2)
Sigma_temp = reader.getIndividualSigma(M=10, d=2)
Vt_temp = reader.getIndividualVt(M=10, d=2)
Pseudoinverse_temp = reader.getIndividualPseudoinverse(M=10, d=2)
#'''


'''
#reads in the full S data matrix
S_data = reader.loadSLookupTables()

#reads a particular S
S_temp = reader.getIndividualS(M=10, d=2)
#'''

W_data = reader.loadWLookupTables()


#gets the W data
W_temp = reader.getIndividualW(M=10, d=2, l=1)

#gets the Y and the Z data
Y_temp, Z_temp = YZReader.loadYZTables()



tomato = 0
