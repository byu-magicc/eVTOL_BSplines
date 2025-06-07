import os, sys
from pathlib import Path


sys.path.insert(0,os.fspath(Path(__file__).parents[1]))
tempPath = sys.path

#imports the lookup table generator
from eVTOL_BSplines.path_generation_helpers.lookup_table_generators import lookUpTablesGenerator


#instantiates the generator for the thing on the other thing
generator = lookUpTablesGenerator()

#calls the function to generate the B lookup tables
generator.generateBLookupTables()

#generates the S Lookup tables
generator.generateSLookupTables()

#generates the W lookup tables
generator.generateWLookupTables()