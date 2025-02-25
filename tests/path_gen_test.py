import os, sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from bsplinegenerator.bsplines import BsplineEvaluation

sys.path.insert(0,os.fspath(Path(__file__).parents[1]))
tempPath = sys.path




potato = 0