import ctypes 
import pathlib 
import os 
import numpy as np
from pathlib import Path

script_dir = os.path.abspath(os.path.dirname(__file__))
libname_str = os.path.join(script_dir)
libname = pathlib.Path(libname_str)



so_file_path = os.path.join(os.fspath(Path(__file__).parents[1]), 
                       "build",
                       "libPathObjectivesAndConstraints.so")
#creates the .so path
lib = ctypes.CDLL(so_file_path)





sweetPotato = 100