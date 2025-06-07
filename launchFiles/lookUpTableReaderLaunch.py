#this is a test file to see if our generator for the lookup tables worked,
#and to see if our readers are working well just the same

import os, sys
from pathlib import Path


sys.path.insert(0,os.fspath(Path(__file__).parents[1]))
tempPath = sys.path

from eVTOL_BSplines.path_generation_helpers.lookup_table_helpers import lookUpTableReader, YZGeneratorReader



#creates the reader
reader = lookUpTableReader()



potato = 0