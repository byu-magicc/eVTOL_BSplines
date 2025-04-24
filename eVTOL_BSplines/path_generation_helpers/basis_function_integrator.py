from matrix_generators_efficient import integrateBasisFunctions
import os, sys
from pathlib import Path



temp1 = os.fspath(Path(__file__).parents[0])
inputFileName = os.path.abspath(os.path.join(temp1, 'lookUpTables/degree_5_basis.npz'))
outputFileName = os.path.abspath(os.path.join(temp1, 'lookUpTables/degree_5_integrations.npz'))

integrator = integrateBasisFunctions(fileName=inputFileName, highestDegree=5)


integrator.createIntegrations(outputFileName=outputFileName)