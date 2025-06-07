import os, sys
from pathlib import Path


sys.path.insert(0,os.fspath(Path(__file__).parents[1]))
tempPath = sys.path

#imports the lookup table generator
from eVTOL_BSplines.path_generation_helpers.lookup_table_helpers import lookUpTablesGenerator, YZGeneratorReader


#instantiates the generator for the thing on the other thing
generator = lookUpTablesGenerator()

#instantiates the YZ generator
YZ_gen = YZGeneratorReader()

#calls 