#this file tests the SFC corridor constraints
import numpy as np
from path_generation.safe_flight_corridor import SFC


#creates the dimension
dimension = np.array([[10.],
                      [20.],
                      [30.]])

#creates the rotation
rotation_CL_to_world = np.array([[1.0, 0.0, 0.0],
                                 [0.0, 1.0, 0.0],
                                 [0.0, 0.0, 1.0]])

rotation_world_to_CL = rotation_CL_to_world.T

#creates the translation vector in the world frame
translation_world = np.array([[0.0],
                              [0.0],
                              [0.0]])

#creates the translatio n vector in the CL frame
translation_CL = rotation_world_to_CL @ translation_world

#creates the SFC
tempSFC = SFC(dimensions=dimension,
              translation=translation_CL,
              rotation=rotation_CL_to_world)

#calls the function to get the normals and vertices
tempSFC.getNormalsVertices_3d()