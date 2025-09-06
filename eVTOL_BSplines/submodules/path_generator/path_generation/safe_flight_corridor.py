import numpy as np
from dataclasses import dataclass
import matplotlib.pyplot as plt

@dataclass
class SFC:
    """Safe Flight Corridor Data Class"""
    # rotation @ translation gives true center of sfc
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
        
        #gets the points and omits the last one (because that's how it works for this thing now)
        vertices = (self.getPointsToPlot())[:,:-1]

        numVertices = np.shape(vertices)[1]

        #creates the matrix of normal vectors
        normalVectors = np.ndarray((2,0))

        #now we iterate over to get the normal vectors 
        for i in range(numVertices):
            #gets the current vertex
            currentVertex = vertices[:,i:(i+1)]

            if i == (numVertices - 1):
                nextVertex = vertices[:,0:1]
            else:
                nextVertex = vertices[:,(i+1):(i+2)]
            

            #gets the vector from the current vertex to the next vertex
            vector_currentToNext = nextVertex - currentVertex

            #gets it normalized
            vector_currentToNext_normal = vector_currentToNext / np.linalg.norm(vector_currentToNext)

            #gets the rotation matrix 
            R_norm = np.array([[0, -1],
                               [1, 0]])
            
            #gets the normal vector
            vector_normal = R_norm @ vector_currentToNext_normal

            #concatenates onto the normalVectors matrix
            normalVectors = np.concatenate((normalVectors, vector_normal), axis=1)


        

        return vertices, normalVectors

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

    #defines the function to get the A Matrix
    

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



