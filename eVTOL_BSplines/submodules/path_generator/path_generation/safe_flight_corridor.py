import numpy as np
from dataclasses import dataclass
import matplotlib.pyplot as plt
from shapely.geometry import MultiPoint


#defines the new list of edge vertices indices
edges_verticesIndices = [[0,1], #pair 0 
                         [0,3], #pair 1 
                         [0,4], #pair 2
                         [1,2], #pair 3
                         [1,5], #pair 4
                         [2,3], #pair 5
                         [2,6], #pair 6
                         [3,7], #pair 7
                         [4,5], #pair 8
                         [4,7], #pair 9
                         [5,6], #pair 10
                         [6,7]] #pair 11


@dataclass
class SFC:
    """Safe Flight Corridor Data Class"""
    # rotation @ translation gives true center of sfc
    #dimensions here does NOT state the number of dimensions (2-D, 3-D, etc)
    #it refers to the size of each of those dimensions for the rectangle or rectangular prism.
    #if it is 2-D, it will be the length and width as an array
    #if it is 3D, it will be the length width and height of the rectangle, etc.
    dimensions: np.ndarray
    translation: np.ndarray
    rotation: np.ndarray

    def getRotatedBounds(self):
        max_bounds = self.translation + self.dimensions/2
        min_bounds = self.translation - self.dimensions/2
        return min_bounds, max_bounds
    
    def getPointsToPlot(self):
        if(len(self.dimensions.flatten()) == 2): #2 D
            return self.getPointsToPlot2D()
        else: # 3 D
            return self.getPointsToPlot3D()
    
    def getPointsToPlot2D(self):
        min_bounds, max_bounds = self.getRotatedBounds()
        x_min = min_bounds.item(0)
        x_max = max_bounds.item(0)
        y_min = min_bounds.item(1)
        y_max = max_bounds.item(1)
        points_unrotated = np.array([[x_min, x_min, x_max, x_max, x_min],
                                     [y_min, y_max, y_max, y_min, y_min]])
        points_rotated = self.rotation @ points_unrotated
        return points_rotated
    
    
    def getPointsToPlot3D(self):
        min_bounds, max_bounds = self.getRotatedBounds()
        x_min = min_bounds.item(0)
        x_max = max_bounds.item(0)
        y_min = min_bounds.item(1)
        y_max = max_bounds.item(1)
        z_min = min_bounds.item(2)
        z_max = max_bounds.item(2)
        points = np.array([[x_max, x_min, x_min, x_max, x_max, x_min, x_min, x_max, x_max, x_max, x_max, x_max, x_min, x_min, x_min, x_min],
                           [y_min, y_min, y_max, y_max, y_max, y_max, y_min, y_min, y_min, y_max, y_max, y_min, y_min, y_min, y_max, y_max],
                           [z_min, z_min, z_min, z_min, z_max, z_max, z_max, z_max, z_min, z_min, z_max, z_max, z_max, z_min, z_min, z_max]])
        points = self.rotation @ points
        return points

    #gets the normal vectors with each corresponding vertex.
    #IMPORTANT NUANCE: Normal vectors are pointing OUTWARD from the convex hull
    #this makes it easier to use the less than or equal operator than 
    def getNormalsVertices_2d(self):
        
        #gets the x dimension
        x_dimension = self.dimensions[0]
        #gets the y dimension
        y_dimension = self.dimensions[1]

        x_min = -x_dimension/2
        x_max = x_dimension/2

        y_min = -y_dimension/2
        y_max = y_dimension/2

        #creates the initial vertices before rotation and translation
        initialVertices = np.array([[x_max, x_max, x_min, x_min],
                                    [y_max, y_min, y_min, y_max]])
        


        #gets the two rotation matrices
        rotation_CL_to_world = self.rotation
        rotation_world_to_CL = rotation_CL_to_world.T

        #gets the rotation of the initial vertices
        rotatedVertices = rotation_CL_to_world @ initialVertices

        #gets the translation
        translation_CL = self.translation
        translation_World = rotation_CL_to_world @ translation_CL

        #adds the tranlsation to the rotated vertices
        self.finalVertices = rotatedVertices + translation_World

        numVertices = np.shape(self.finalVertices)[1]


        #creates the 90 degree rotation matrix to the right
        Rotation_Right = np.array([[0.0, -1.0],
                                   [1.0, 0.0]])


        normalVectors_list = []

        vertices_list = []

        for i in range(numVertices):

            currentVertex = self.finalVertices[:,i:(i+1)]

            #appends the current vertex to the vertices list
            vertices_list.append(currentVertex)

            if i == (numVertices - 1):
                nextVertex = self.finalVertices[:,0:1]
            else:
                nextVertex = self.finalVertices[:,(i+1):(i+2)]

            #gets the vector from current to next
            vectorCurrentToNext = nextVertex - currentVertex

            #gets the normal of that
            vectorCurrentToNext_norm = vectorCurrentToNext / np.linalg.norm(vectorCurrentToNext)

            #gets the normal vector
            currentNormalVector = Rotation_Right @ vectorCurrentToNext_norm

            normalVectors_list.append(currentNormalVector)

    
        return normalVectors_list, vertices_list


    #does the same thing for the 3-Dimensional case
    def getNormalsVertices_3d(self):
        #gets the size of the x dimension
        x_dimension = self.dimensions[0,0]
        #gets the size of the y dimensions
        y_dimension = self.dimensions[1,0]
        #same for z
        z_dimension = self.dimensions[2,0]

        #gets the mins and maxes for the x y and z
        x_min = -x_dimension/2.0
        x_max = x_dimension/2.0

        y_min = -y_dimension/2.0
        y_max = y_dimension/2.0

        z_min = -z_dimension/2.0
        z_max = z_dimension/2.0

        initialVertices = np.array([[x_min, x_max, x_max, x_min, x_min, x_max, x_max, x_min],
                                         [y_min, y_min, y_max, y_max, y_min, y_min, y_max, y_max],
                                         [z_min, z_min, z_min, z_min, z_max, z_max, z_max, z_max]])

        #gets the two rotation matrices
        rotation_CL_to_world = self.rotation
        rotation_world_to_CL = rotation_CL_to_world.T

        #gets the rotated vertices
        rotatedVertices = rotation_CL_to_world @ initialVertices

        #gets the translation in the CL frame
        translation_CL = self.translation
        #gets the translatio nin the world grame
        translation_World = rotation_CL_to_world @ translation_CL

        #adds the translation t othe rotated vertices
        self.finalVertices = rotatedVertices + translation_World

        #calls the function to get the 3d vertices and nromals
        self.normalVectors_list, self.verticesForNormalVectors_list = self.generate3DNormalsVertices(vertices=self.finalVertices)

        return self.normalVectors_list, self.verticesForNormalVectors_list

    #defines the function to get all vertices
    def getAllVertices(self):
        self.getNormalsVertices_3d()
        return self.finalVertices

    def getAllVerticesList(self):
        
        verticesShape = np.shape(self.finalVertices)
        numVertices = verticesShape[1]

        allVerticesList = []

        for i in range(numVertices):
            currentVertex = (self.finalVertices)[:,i:(i+1)]
            allVerticesList.append(currentVertex)

        return allVerticesList
    
    #function to get the A and b matrices
    def getAbMatrices(self):

        normals, vertices = self.getNormalsVertices()

        numDimensions = np.size(self.dimensions)

        self.A = np.ndarray((0, numDimensions))

        self.b = np.ndarray((0, 1))

        for normal, vertex in zip(normals, vertices):

            self.A = np.concatenate((self.A, normal.T), axis=0)
            
            b_value = normal.T @ vertex
            self.b = np.concatenate((self.b, b_value), axis=0)

        return self.A, self.b

    #gets the normals and vertices
    def getNormalsVertices(self):
        #checks the dimensionality of the problem
        numDimensions = np.size(self.dimensions)
        
        #if the number of dimensions is 2 or 3, we get the 2d or 3d normals and vertices
        if numDimensions == 2:
            return self.getNormalsVertices_2d()
        elif numDimensions == 3:
            return self.getNormalsVertices_3d()
        


    #defines the function to get the shapely convex hull
    def getConvexHull(self):

        #gets the vertices list
        verticesList = self.getAllVerticesList()
        #constructs the convex hull
        convex_hull = MultiPoint(verticesList).convex_hull

        return convex_hull


    #generates normals and vertices for the 3d case

    #defines the function to get the rotation matrix
    #which rotates from corridor frame to world frame
    def getRotation_corridorToWorld(self):
        return self.rotation
    
    #defines the function to get the rotation matrix
    #which rotates from the world frame to the corridor frame
    def getRotation_worldToCorridor(self):
        return self.rotation.T


    #defines a function to set the old vertices and normal vectors
    def setNormalsVertices_old(self,
                               vertices: list[np.ndarray],
                               normalVectors: list[np.ndarray]):
        
        self.vertices_old = vertices
        self.normalVectors_old = normalVectors
    
    def getNormalsVertices_old(self):
        return self.normalVectors_old, self.vertices_old

    #defines the function to get the normal vectors in 3D
    def generate3DNormalsVertices(self,
                                  vertices: np.ndarray):


        #creates the list of the edge vectors
        self.edgeLengths = []
        self.edgeVectors = []
        #creates the end
        for edgeVertices in edges_verticesIndices:

            #gets the indices for the two vertices
            startVertexIndex = edgeVertices[0]
            endVertexIndex = edgeVertices[1]

            #gets the corresponding vertices
            startVertex_pos = vertices[:,startVertexIndex:(startVertexIndex+1)]
            endVertex_pos = vertices[:,endVertexIndex:(endVertexIndex+1)]

            #gets the edge vector
            currentEdgeVector = endVertex_pos - startVertex_pos
            
            currentEdgeVector_length = np.linalg.norm(currentEdgeVector)
            #gets it normalized
            currentEdgeVector_norm = currentEdgeVector / currentEdgeVector_length
            #appends the unit vector
            self.edgeVectors.append(currentEdgeVector_norm)
            self.edgeLengths.append(currentEdgeVector_length)


        edgeVectorPairsIndices_list = [[1,0],
                                       [0,2],
                                       [3,4],
                                       [5,6],
                                       [2,1],
                                       [8,9]]                                

        #creates the normals vectors list
        normalVectorsList = []

        #creates the list of the points that correspoints to each normal vector
        verticesForNormalVectorsList = []

        #ok. now with the edge vectors list, we can go back through and use cross products 
        # that define the exterior-facing vectors
        for edgeVectorPairIndices in edgeVectorPairsIndices_list:

            #gets the corresponding position, which is just the position
            #of the 

            #gets the primary vector
            primaryVector_index = edgeVectorPairIndices[0]
            secondaryVector_index = edgeVectorPairIndices[1]

            #Gets the primary and the secondary vectors
            primaryVector = self.edgeVectors[primaryVector_index]
            secondaryVector = self.edgeVectors[secondaryVector_index]

            #gets the primary vector shape
            primaryVector_shape = np.shape(primaryVector)

            #gets the flattened versions of the primary and secondary vectors
            primaryVector_flattened = primaryVector.flatten()
            secondaryVector_flattened = secondaryVector.flatten()

            #gets the cross product of the primary and secondary vectors
            normalVector_temp_flattened = np.cross(primaryVector_flattened, secondaryVector_flattened)
            #reshapes the normal vector
            normalVector_temp = np.reshape(normalVector_temp_flattened, primaryVector_shape)
            #appends to the normal Vectors list
            normalVectorsList.append(normalVector_temp)

            #gets the edge vertices for the primary vector
            primaryVectorVertices_indices = edges_verticesIndices[primaryVector_index]
            #gets the start position for the primary vector
            primaryVector_startPositionIndex = primaryVectorVertices_indices[0]
            #finally gets the primary vector start position
            primaryVector_startPosition = vertices[:,primaryVector_startPositionIndex:(primaryVector_startPositionIndex+1)]
            #appends to the vertices for normal vectors list
            verticesForNormalVectorsList.append(primaryVector_startPosition)


        return normalVectorsList, verticesForNormalVectorsList
    

    def getEdgeVectors(self):
        self.getNormalsVertices()
        return self.edgeVectors, self.edgeLengths, edges_verticesIndices
    

class SFC_Data:
    def __init__(self, sfc_list: list, point_sequence: np.ndarray, min_num_intervals_per_corridor: int = 1):
        self._sfc_list = sfc_list
        # print("sfc list: ", self._sfc_list)
        self._num_corridors = len(self._sfc_list)
        # print("num corridors: " , self._num_corridors)
        self._point_sequence = point_sequence
        self._min_num_intervals_per_corridor = min_num_intervals_per_corridor
        self._intervals_per_corridor = self.__evaluate_intervals_per_corridor()
        self._num_intervals = np.sum(self._intervals_per_corridor)

    def get_sfc_list(self)->list[SFC]:
        return self._sfc_list

    def get_num_corridors(self):
        return self._num_corridors
    
    def get_point_sequence(self):
        return self._point_sequence
    
    def get_intervals_per_corridor(self):
        return self._intervals_per_corridor
    
    def get_num_intervals(self):
        return self._num_intervals
    
    def __evaluate_intervals_per_corridor(self):
        if self._num_corridors < 2:
            intervals_per_corridor = 5
        else:
            distances = np.linalg.norm(self._point_sequence[:,1:] - self._point_sequence[:,0:-1],2,0)
            min_distance = np.min(distances)
            intervals_per_corridor = []
            for i in range(self._num_corridors):
                num_intervals = (int(np.round(distances[i]/min_distance)))*self._min_num_intervals_per_corridor
                intervals_per_corridor.append(num_intervals)
        # intervals_per_corridor = [1,2,3]
        # intervals_per_corridor = [2,2,3]
        # intervals_per_corridor = [2,3,4]
        # intervals_per_corridor = [2,2,2,2,6]
        return intervals_per_corridor

def plot_sfc(sfc: SFC, ax,alpha=1):
    if(len(sfc.dimensions.flatten()) == 2): #2D
        plot_2D_sfc(sfc, ax, alpha)
    else: # 3D
        plot_3D_sfc(sfc, ax, alpha)

def plot_sfcs(sfcs: list, ax, alpha=1):
    if sfcs != None:
        for sfc_index in range(len(sfcs)):
            plot_sfc(sfcs[sfc_index], ax, alpha)

def plot_2D_sfc(sfc: SFC, ax, alpha=1):
    points_rotated = sfc.getPointsToPlot()
    ax.plot(points_rotated[0,:], points_rotated[1,:], alpha=alpha)

def plot_3D_sfc(sfc: SFC, ax, alpha=1):
    points_rotated = sfc.getPointsToPlot()
    ax.plot(points_rotated[0,:], points_rotated[1,:],points_rotated[2,:], alpha= alpha)

def get2DRotationAndTranslationFromPoints(point_1,point_2):
    # returns rotation transforms x_vector to vector paralell
    distance = point_2 - point_1
    dx = distance.item(0)
    dy = distance.item(1)
    psi = np.arctan2(dy,dx)
    c_psi = np.cos(psi)
    s_psi = np.sin(psi)
    rotation = np.array([[c_psi, -s_psi],
                 [s_psi, c_psi]])
    translation = rotation.T @ (point_1 + point_2)/2
    min_length = np.linalg.norm(distance,2)
    return rotation, translation, min_length

def get3DRotationAndTranslationFromPoints(point_1,point_2):
    # returns rotation transforms x_vector to vector paralell
    distance_1 = point_2 - point_1
    dx_1 = distance_1.item(0)
    dz_1 = distance_1.item(2)
    theta = np.arctan2(dz_1,dx_1)
    c_theta = np.cos(theta)
    s_theta = np.sin(theta)
    Ry = np.array([[c_theta, 0, s_theta],
                [0, 1, 0],
                [-s_theta, 0, c_theta]])
    distance_2 = Ry @ distance_1
    dx_2 = distance_2.item(0)
    dy_2 = distance_2.item(1)
    psi = np.arctan2(dy_2,dx_2)
    c_psi = np.cos(psi)
    s_psi = np.sin(psi)
    Rz = np.array([[c_psi, -s_psi, 0],
                 [s_psi, c_psi, 0],
                 [0,       0,       1]])
    rotation = Ry.T @ Rz
    translation = rotation.T @ (point_1 + point_2)/2
    min_length = min_length = np.linalg.norm(distance_1,2)
    return rotation, translation, min_length




