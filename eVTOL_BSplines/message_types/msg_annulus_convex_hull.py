import numpy as np

#this file implements the class to communicate the convex hulls for a circle.
#or rather the partial negative space of a circle, which is really called an annulus
#and the problem with an annulus is that it is not a convex hull.
#But........ If turn it into a bunch of pizza slices, and cut off the round parts, then it is a convex hull

#sets the north reference vector
northReference = np.array([[1.0],[0.0]])

class Msg_Annulus_Convex_Hull:

    #creates the initialization function
    #arguments:
    #centerPosition: the point of the center of the circle
    #innerRadius: the smaller radius for the annulus
    #outerRadius: the outer radius for the annulus
    #centerToPivot_unit: the unit vector that points from the center to the pivot point
    #arcAngle: the angle over which there are the valid convex hulls. centerAngle position bisects this
    def __init__(self,
                 centerPosition: np.ndarray,
                 innerRadius: float,
                 outerRadius: float,
                 centerToPivot_unit: np.ndarray,
                 arcAngle: float,
                 numSections: int):
        
        self.northUnit = np.array([[1.0],
                              [0.0]])

        #saves the center position here
        self.centerPosition = centerPosition
        
        self.innerRadius = innerRadius
        self.outerRadius = outerRadius

        self.centerToPivot_unit = centerToPivot_unit

        #gets the center Angle
        self.centerAngle = getAngleVectors(v1=self.northUnit,
                                           v2=self.centerToPivot_unit)

        self.arcAngle = arcAngle

        self.numSections = numSections

        #calls the helper function to get the angles list
        self.getAnglesLists()

        #calls the helper function to get the hull vertices list
        self.getHullVerticesNormals()
        

    #defines the function to get the angles lists
    def getAnglesLists(self):
        #gets the start and end angles, moving clockwise around the circle
        self.startAngle = self.centerAngle - self.arcAngle / 2.0
        self.endAngle = self.centerAngle + self.arcAngle / 2.0

        #gets the individual angle of each bin
        self.incrementalAngle = self.arcAngle / self.numSections

        #creates the list of each angle for each bin (num angles = num bins + 1)
        anglesList = [self.startAngle]

        tempAngle = self.startAngle
        #iterates over tall the angles
        for i in range(self.numSections):
            tempAngle += self.incrementalAngle
            anglesList.append(tempAngle)

        #innerradius list
        self.innerRadiusPosition_list = []
        self.outerRadiusPosition_list = []

        for angle in anglesList:
            
            #gets the R
            R_temp = getRotationMatrix(theta=angle)
            #gets the innerRadius position
            innerRadius_position = self.centerPosition + self.innerRadius * R_temp @ self.northUnit
            outerRadius_position = self.centerPosition + self.outerRadius * R_temp @ self.northUnit

            self.innerRadiusPosition_list.append(innerRadius_position)
            self.outerRadiusPosition_list.append(outerRadius_position)


    #and now with the list of angles, we can construct the list of coordinates
    #as in obtain the list of vertices
    def getHullVerticesNormals(self):

        self.allVerticesList = []
        self.allNormalsList = []

        for i in range(len(self.innerRadiusPosition_list) - 1):

            #gets the inner position 1
            innerPosition_1 = self.innerRadiusPosition_list[i]
            innerPosition_2 = self.innerRadiusPosition_list[i+1]
            outerPosition_1 = self.outerRadiusPosition_list[i]
            outerPosition_2 = self.outerRadiusPosition_list[i+1]

            #puts together the vertices list, in the following order: outerRight, outerLeft, innerLeft, innerRight
            verticesList = [outerPosition_2, outerPosition_1, innerPosition_1, innerPosition_2]

            normalsList = getNormalVectors(verticesList=verticesList)
            
            #appends to the vertices and the normals
            self.allVerticesList.append(verticesList)
            self.allNormalsList.append(normalsList)

    #creates the function to return the list of hull vertices
    def getVertices(self):
        return self.allVerticesList
    
    def getNormals(self):
        return self.allNormalsList


        
def getRotationMatrix(theta: float):

    R = np.array([[np.cos(theta), -np.sin(theta)],
                  [np.sin(theta), np.cos(theta)]])
    
    return R
    

#defines the function to get the angle between two vectors
def getAngleVectors(v1: np.ndarray,
                    v2: np.ndarray)->float:
    
    #gets  the dimension of the vector
    if np.size(v1) == 3:
        v1 = v1[:-1,0:1]
        v2 = v2[:-1,0:1]

    #gets the cross product
    normals_cross = v1[0,0]*v2[1,0] - v1[1,0]*v2[0,0]

    #gets the dot product
    normals_dot = v1[0,0]*v2[0,0] + v1[1,0]*v2[1,0]

    #gets the theta
    theta = np.arctan2(normals_cross, normals_dot)

    return theta


#defines the function to obtain the normal vectors from the vertices list
def getNormalVectors(verticesList: list[np.ndarray])->list[np.ndarray]:
    
    numVertices = len(verticesList)

    normalVectors = []

    for i in range(numVertices):

        #gets the start vertex and the end vertex
        startVertex = verticesList[i]
        #gets the end vertex
        endVertex = verticesList[(i+1)%numVertices]

        #gets the current edge
        currentEdge = endVertex - startVertex
        #normalizes it
        currentEdge_norm = currentEdge / np.linalg.norm(currentEdge)
        #gets the rotation matrix
        R = getRotationMatrix(theta=np.pi/2)
        #gets the current normal vector
        currentNormal = R @ currentEdge_norm
        normalVectors.append(currentNormal)

    #returns the normal vectors list
    return normalVectors
