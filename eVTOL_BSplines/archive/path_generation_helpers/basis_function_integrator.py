from matrix_generators_efficient import integrateBasisFunctionsContinuous
import os, sys
from pathlib import Path
import numpy as np


#creates the temp file path for the degree 5 integrations
filePath = 'lookUpTables2/degree_5_integrations.npz'



temp1 = os.fspath(Path(__file__).parents[0])
outputFileName = os.path.abspath(os.path.join(temp1, filePath))

integrator = integrateBasisFunctionsContinuous(outputFileName=outputFileName, highestDegree=5)





loadedDictionary = np.load(outputFileName)
list_of_arrays = [loadedDictionary[f'degree_{i}'] for i in range(len(loadedDictionary.files))]



tamale = 0