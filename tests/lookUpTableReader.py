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

d = 2

M = 5

ell = d - 1

key = "degree2_l1_M5"

degree, ell, M = reader.getIndividualWMetadata(key=key)


'''
W = reader.getIndividualW(d=d, ell=ell, M=M)

S = reader.getIndividualS(d=d, M=M)

B, U1, U2, Sigma, Vt = reader.getIndividualB(d=d, M=M)

#'''

potato = 0