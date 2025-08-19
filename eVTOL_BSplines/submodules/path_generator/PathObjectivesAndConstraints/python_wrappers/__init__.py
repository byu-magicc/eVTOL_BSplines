import os
import ctypes
from pathlib import Path



lib_path = os.path.join(
    os.fspath(Path(__file__).parents[1]), 
    "build", 
    "libPathObjectivesAndConstraints.so"
)

# Load the shared library
constraints_lib = ctypes.CDLL(lib_path)

from .curvature_constraints import *
from .incline_constraints import *
from .objective_functions import *
from .obstacle_constraints import *
from .waypoint_constraints_old import *
from .waypoint_constraints import *
