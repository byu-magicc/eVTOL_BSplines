import numpy as np
from path_generation.safe_flight_corridor import SFC_Data, SFC
import cvxpy as cp
from bsplinegenerator.bspline_to_minvo import convert_list_to_minvo_control_points



#defines the function to get the number of control points from an existing array
def getNumCtrPts_array(controlPoints: np.ndarray):

    #gets the shape of the control points array
    controlPoints_shape = np.shape(controlPoints)

    #case it is along the zero axis
    if controlPoints_shape[0] > controlPoints_shape[1]:
        numControlPoints = controlPoints_shape[0]
    else:
        numControlPoints = controlPoints_shape[1]

    #returns the number of control points
    return numControlPoints


#gets the central control points from the complete initial control points list
def getPartitionedControlPoints_initial(initialControlPoints: np.ndarray,
                                        degree: int):
    
    #checks the shape and works with that shape
    pointsShape = np.shape(initialControlPoints)

    #case tall matrix
    if pointsShape[0] > pointsShape[1]:
        numControlPoints = pointsShape[0]
        numIntervalsOfInterest = numControlPoints - degree
        startControlPoints_initial = initialControlPoints[0:degree,:]
        centralControlPoints_initial = initialControlPoints[degree:numIntervalsOfInterest,:]
        endControlPoints_initial = initialControlPoints[numIntervalsOfInterest:,:]
    #case fat matrix
    else:
        numControlPoints = pointsShape[1]
        numIntervalsOfInterest = numControlPoints - degree
        startControlPoints_initial = initialControlPoints[:,0:degree]
        centralControlPoints_initial = initialControlPoints[:,degree:numIntervalsOfInterest]
        endControlPoints_initial = initialControlPoints[:,numIntervalsOfInterest:]

    #returns the initial central control points
    return startControlPoints_initial, centralControlPoints_initial, endControlPoints_initial

#function to reconstruct flattened control points, with the start and end sections
def reconstructFlattenedControlPoints(startControlPoints: np.ndarray,
                                      flattenedCenterControlPoints: np.ndarray,
                                      endControlPoints: np.ndarray,
                                      numCenterControlPoints: int,
                                      dimension: int):
    
    #reshapes the flattened center control points into the correct shape
    centerControlPoints = np.reshape(flattenedCenterControlPoints, (dimension, numCenterControlPoints))

    #concatenates them together to create the full control point array
    completeControlPoints = np.concatenate((startControlPoints, centerControlPoints, endControlPoints), axis=1)

    #returns  the control points
    return completeControlPoints

#defines a function to convert from an sfc list to cvxpy bounds
#Arguments: 
#sfc_list: the list of safe flight corridors for this thing
#controlPoints_parsed: control points parsed so that they represent the correct numbers
def sfcListToCvxpyConstraints(sfc_list: list[SFC],
                              controlPoints_parsedLengths: list[int],
                              cpPoints: cp.Variable):
    
    #creates the cpPoint index
    cpPoint_index = 0
    #creates the constraints array
    constraints = []
    #iterates over each sfc list
    for sfc, cntPts_parsedLength in zip(sfc_list, controlPoints_parsedLengths):
        
        #gets the temp rotation
        tempRotation_CorridorToWorld = sfc.rotation
        #gets the rotation from World to corridor frame
        tempRotation_WorldToCorridor = tempRotation_CorridorToWorld.T
        #gets the translation
        tempTranslation = sfc.translation
        #gets the x dimension
        dimensions = sfc.dimensions
        x_dimension = dimensions[0]
        y_dimension = dimensions[1]

        #gets the locale, which is 


        #iterates over all the cp points belonding to this partcular
        for i in range(cntPts_parsedLength):
            tempIndex = cpPoint_index + i
            #gets the current point
            tempCpPoint = cpPoints[:, tempIndex]
            #now, let's set the boundaries
            locale = tempRotation_WorldToCorridor @ (tempCpPoint - tempTranslation)

            #sets the boundaries for the locale for the temp constraints
            tempConstraints = \
            [locale[0] <= x_dimension/2.0,
            locale[0] >= -x_dimension/2.0,
            locale[1] <= y_dimension/2.0,
            locale[1] >= -y_dimension/2.0]

            #adds the temp constraints to the full list of constraints
            constraints += tempConstraints


        #adds the parsed lengths to the cp index
        cpPoint_index += cntPts_parsedLength

    #returns the constraints
    return constraints


#defines a function to get a list of the number of minvo points associated with a list of certain minvo points
def getNumMinvoPoints_parsed(controlPointsParsedNumberList: list[int],
                             degree: int):

    numMinvoPointsList = []
    #iterates over each number
    for numControlPoints in controlPointsParsedNumberList:

        #gets the number of minvo points on a section
        numMinvoSections = numControlPoints - degree
        numMinvoPoints = numMinvoSections * (degree+1)

        #saves this to the num minvo points list
        numMinvoPointsList.append(numMinvoPoints)

    #returns the num minvo points list
    return numMinvoPointsList


#defines the function to get the minvo points 
def controlPointParsed_toMinvo(controlPointParsedList: list[np.ndarray],
                               degree: int):


    minvoPointsList = []
    #iterates over the control points list
    for controlPointsSection in controlPointParsedList:

        #gets the minvo points section
        minvo_points_section = convert_list_to_minvo_control_points(bspline_control_points=controlPointsSection,
                                                                      order=degree)


        minvoPointsList.append(minvo_points_section)

    return minvoPointsList


#defines a function to generate the b matrix 
def generate_A_b(normalVectors: list[np.ndarray],
                 vertices: list[np.ndarray])->tuple[np.ndarray, np.ndarray]:
    
    #this function assumes that the shape of the normal vectors and the vertices arrays
    #are given as (dimension x N)

    normalVectors = np.concatenate(normalVectors, axis=1)
    vertices = np.concatenate(vertices, axis=1)

    A = normalVectors.T

    #gets the number of vectors here
    numVectors = np.shape(A)[0]

    #creates the b vector
    b = np.zeros((numVectors,1))

    for i in range(numVectors):

        #gets the n transposed vector
        n_transpose_temp = A[i:(i+1),:]

        #gets the vertex
        vertex_temp = vertices[:,i:(i+1)]

        #gets the temp result
        tempResult = n_transpose_temp @ vertex_temp

        #and then indexes it
        tempResult = tempResult[0,0]

        #and then saves it
        b[i,0] = tempResult

    #returns the A and b matrix
    return A, b


#defines the function to generate the list of A and b matrices from the SFC list
def generateConstraintsSFC(sfc_data: SFC_Data):

    #gets the list of sfcs
    sfc_list = sfc_data.get_sfc_list()

    #creates the list for the A and b matrices
    A_list = []

    b_list = []

    #iterates over each sfc in the list
    for sfc_temp in sfc_list:

        #gets the normal vectors and vertices
        normalVectors_temp, vertices_temp = sfc_temp.getNormalsVertices_2d()

        #gets the A and b matrices
        A_temp, b_temp = generate_A_b(normalVectors=normalVectors_temp,
                                      vertices=vertices_temp)
        
        A_list.append(A_temp)

        b_list.append(b_temp)

    #returns the A and b lists
    return A_list, b_list


#defines the function to get the shortened control points
def getShortenedControlPoints(startControlPoints: np.ndarray,
                              endControlPoints: np.ndarray):
    
    outputStartControlPoints = startControlPoints[:2,:]
    outputEndControlPoints = endControlPoints[:2,:]

    return outputStartControlPoints, outputEndControlPoints


#defines the function to get the list of number of control points
#for each respective section of the whole thing
def getNumCntPts_list(sfc_data: SFC_Data,
                      numPointsPerUnit: float)->list[int]:

    #gets the sfc_list
    sfc_list = sfc_data.get_sfc_list()

    numCntPts_list = []

    for sfc in sfc_list:


        #gets the x dimension
        x_dimension = sfc.dimensions[0]

        #gets the estimated number of control points
        numControlPoints = int(x_dimension*numPointsPerUnit) + 1

        numCntPts_list.append(numControlPoints)
        
    #returns the list of the number of control points
    return numCntPts_list


#defines the function to get the total number of control points
def getNumCntPts(numCntPts_list: list[int],
                degree: int):
    #get the initial sum
    initialSum = sum(numCntPts_list)
    #gets the number of corridors, which is the number of items in the list
    numCorridors = len(numCntPts_list)
    #gets the num control points
    numCntPts = initialSum - degree * (numCorridors - 1)
    #returns the total number of contorl points
    return numCntPts


#stacks the columns of a matrix up
def stackMatrix(controlPoints: np.ndarray):

    #gets the shape of control points
    controlPoints_shape = np.shape(controlPoints)
    
    controlPoints_stacked = np.ndarray((0,1))
    for i in range(controlPoints_shape[1]):
        #gets the control points partition
        controlPointsPartition = controlPoints[:,i:(i+1)]

        #stacks it up
        controlPoints_stacked \
            = np.concatenate((controlPoints_stacked, controlPointsPartition), axis=0)
        

    #returns the contorl points stackes
    return controlPoints_stacked




#defines the function to create the arc for two sfcs
def generateArc(sfc_1: SFC,
                sfc_2: SFC):
    



    potato = 0
