
import subprocess
import os
import shutil


#creates the build_directory
path_constraints_directory = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "eVTOL_BSplines",
    "submodules",
    "path_generator",
    "PathObjectivesAndConstraints"
)


potato = path_constraints_directory