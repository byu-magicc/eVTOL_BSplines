#creates a message class to define a safe flight corridor
import numpy as np
from dataclasses import dataclass

#creates the Rotation vector for the 2d normal
R_2d_norm = np.array([[0,-1],
                      [1, 0]])

#sets the number of vertices per face
numVerticesPerFace = 4

class Msg_SFC:
    def __init__(self,
                 dimensions: np.ndarray,
                 translation: np.ndarray,
                 rotation: np.ndarray):
        #create the dimensions of the of the SFC, which is a:
        #(numDimensions x 1) numpy vector. 
        self.dimensions = dimensions

        #translation: the vector from the origin to the center point of the SFC
        #(numDimensions x 1)
        self.translation = translation
        #rotation: the rotation matrix to rotate the sfc to the desired position in space
        #(numDimensions x numDimensions)
        self.rotation = rotation
        
        #sets the dimensionality of the thing
        if len(self.dimensions.flatten()) == 2:
            self.numDimensions = 2
        elif len(self.dimensions.flatten()) == 3:
            self.numDimensions = 3

        #calls the generate characteristics function
        self.generateCharacteristics()

    #defines the function to generate everything for the 2d case
    def generateCharacteristics(self):
        #case we have 2 dimensions in this thing
        if self.numDimensions == 2:
            #vertices section
            #####################################################################
            #first section, we need to creates the rectangle,
            x_max = self.dimensions.item(0) / 2.0
            x_min = -x_max
            y_max = self.dimensions.item(1) / 2.0
            y_min = -y_max

            #creates the initial, unrotated and untranslated vertices
            initialVertices = [np.array([[x_max],[y_max]]), 
                               np.array([[x_max],[y_min]]), 
                               np.array([[x_min],[y_min]]), 
                               np.array([[x_min],[y_max]])]

            #gets the rotated vertices
            rotatedVertices = [self.rotation @ v for v in initialVertices]

            #gets final the translated and rotated vertices
            self.vertices = [self.translation + v for v in rotatedVertices]
            #####################################################################


            #####################################################################
            #edges section to get the edges of the rectangle
            numVertices = len(self.vertices)

            #creates the list of the edges
            self.edges = []

            for i in range(numVertices):
                #gets the current index
                index = i % numVertices 
                #gets the current vertex
                currentVertex = (self.vertices)[index]
                #gets the next index    
                nextVertex = (self.vertices)[index + 1]
                #gets the edge
                currentEdge = nextVertex - currentVertex
                self.edges.append(currentEdge)
            #####################################################################


            #####################################################################
            #normal vectors section
            numEdges = len(self.edges)

            #creates the list of the normals
            self.normals = []

            for i in range(numEdges):
                #gets the temp rotated vector
                tempNorm = R_2d_norm @ (self.edges)[i]

                #appends this to the list
                self.normals.append(tempNorm)

            #####################################################################

        #case we have 3 dimensions in this thing
        elif self.numDimensions == 3:

            #vertices section
            #####################################################################
            x_max = self.dimensions.item(0) / 2.0
            x_min = -x_max
            y_max = self.dimensions.item(1) / 2.0
            y_min = -y_max
            z_max = self.dimensions.item(2) / 2.0
            z_min = -z_max

            #creates the vertices of the 3d SFC. The order Matters, so I'm going to lay this out here
            #if X is Up, Y is to the right, and Z is into the page (the positive axes that is),
            #We start with the coorinate that's completely positive, so thats the top, right, into the page one
            #then we go counterclockwise, at a constant Z, for vertices numbers 0 to 3. Then we repeat this
            #for the face that's in the negative z region
            initialVertices = [np.array([[x_max],[y_max],[z_max]]),#index 0,
                               np.array([[x_max],[y_min],[z_max]]),#index 1,
                               np.array([[x_min],[y_min],[z_max]]),#index 2,
                               np.array([[x_min],[y_max],[z_max]]),#index 3,
                               np.array([[x_max],[y_max],[z_min]]),#index 4,
                               np.array([[x_max],[y_min],[z_min]]),#index 5,
                               np.array([[x_min],[y_min],[z_min]]),#index 6,
                               np.array([[x_min],[y_max],[z_min]])]#index 7,
            
            #gets the rotated vertices
            rotatedVertices = [self.rotation @ v for v in initialVertices]

            #gets the translated vertices
            self.vertices = [self.translation + v for v in rotatedVertices]
            #####################################################################



            #edges section
            #####################################################################
            #for the edges, we go counterclockwise for each of the faces, and I will define the
            #order for each face in the following list
            verticiesListOrder = [[0,3,2,1],#front face
                                  [0,4,7,3],
                                  [0,1,5,4],
                                  [1,2,6,5],
                                  [2,3,7,6],
                                  [4,5,6,7]]
            

            #iterates over the main list
            for verticesSublist in verticiesListOrder:

                #iterates over the verticesSublist
                for i, vertexIndex in enumerate(verticesSublist):
                    
                    next_i = (i+1) % numVerticesPerFace

                    
                    
                    



            #####################################################################


    #creates the function to get the vertices for the sfc
    def get_vertices(self):

        return self.vertices




#The above list was for just one safe flight corridor. 
#next, we need to create a list of safe flight corridors to traverse
class Msg_SFC_list:

    def __init__(self):

        self.corridorList: list[Msg_SFC] = []
        
        pass

    #defines the function to add a corridor to the list
    def addCorridor(self,
                    newCorridor: Msg_SFC):
        
        self.corridorList.append(newCorridor)

    #defines the function to get a corridor list
    def getCorridorList(self):

        return self.corridorList