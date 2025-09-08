#this file re implements the Path Generator's stuff, but this time rewritten to make more sense
#and to be a greater deal of optimality.


import os, sys
import numpy as np
from scipy.optimize import minimize, Bounds, LinearConstraint, NonlinearConstraint, Bounds
from path_generation.matrix_evaluation import get_M_matrix, evaluate_point_on_interval
from PathObjectivesAndConstraints.python_wrappers.objective_functions import ObjectiveFunctions
from PathObjectivesAndConstraints.python_wrappers.curvature_constraints import CurvatureConstraints
from PathObjectivesAndConstraints.python_wrappers.obstacle_constraints import ObstacleConstraints
from PathObjectivesAndConstraints.python_wrappers.incline_constraints import InclineConstraints
from PathObjectivesAndConstraints.python_wrappers.waypoint_constraints import WaypointConstraints
from bsplinegenerator.bspline_to_minvo import get_composite_bspline_to_minvo_conversion_matrix
from path_generation.safe_flight_corridor import SFC_Data, SFC
from path_generation.obstacle import Obstacle
from path_generation.waypoint_data import Waypoint, WaypointData
from eVTOL_BSplines.message_types.msg_control_points import MSG_Control_Points
import time
from bsplinegenerator.bspline_to_minvo import convert_list_to_minvo_control_points


import cvxpy as cp




class SFC_PathGenerator:


    #creates the initialization function.
    def __init__(self,
                 num_intervals_free_space: int = 5):
        #sets the dimension to two for now
        self._dimension = 2
        #sets the number of intervals free space
        self.num_intervals_free_space = num_intervals_free_space
        #I'm only dealing with 3rd degree bsplines (called order in his stuff. I know - it's annoying)
        self._order = 3
        #gets the Minvo conversion matrix
        self._Minvo_matrix = get_M_matrix(order=self._order)
        #gets the objective functions from the c files
        self._objective_func_obj = ObjectiveFunctions(self._dimension)
        self._curvature_const_obj = CurvatureConstraints(self._dimension)
        self._waypoint_const_obj = WaypointConstraints(self._dimension)
        self._obstacle_cons_obj = ObstacleConstraints(self._dimension)



    #defines the function to generate a path through a list of SFCs with direct control point method.
    #that is, we directly modify and optimize the control points and see how that affects things
    #TODO
    def generatePath__directCtrlPts(self,
                                    numIntervalsOfInterestPerCorridor: int,
                                    initialControlPoints: MSG_Control_Points,
                                    sfc_data: SFC_Data = None,
                                    objective_function_type: str = "minimal_velocity_path"):
        
        #gets the initial control points: both the stitched together ones, and the individual ones
        controlPoints_parsed = initialControlPoints.getControlPointsArray_parsed()
        controlPoints_whole = initialControlPoints.getControlPointsArray_complete()
        controlPoints_list = initialControlPoints.getControlPointsArray_list()
        controlPoints_parsedLengths = initialControlPoints.getParsedLengths()

        #gets the initial control points for the whole path and the final ones (3 in each category)
        #these ones will be staying constant the entire time through.
        wholePath_startControlPoints = controlPoints_parsed[0]
        wholePath_endControlPoints = controlPoints_parsed[-1]

        #IMPORTANT. Shaves off altitude. UGGGGHHHHHH!. I keep shooting myself in the foot
        wholePath_startControlPoints = wholePath_startControlPoints[:-1,:]
        wholePath_endControlPoints = wholePath_endControlPoints[:-1,:]

        #gets the control points parsed in the center
        controlPoints_parsed_center = controlPoints_parsed[1:-1]


        #gets the minvo points center and parsed
        minvoPonts_parsed_center = controlPointParsed_toMinvo(controlPointParsedList=controlPoints_parsed_center,
                                                              degree=self._order)

        #gets the applicable parsed lengths matrix (that is one excluding the lengths of the ends)
        controlPoints_parsedLengths_center = controlPoints_parsedLengths[1:-1]

        minvoPoints_parsedLengths_center = getNumMinvoPoints_parsed(controlPointsParsedNumberList=controlPoints_parsedLengths_center,
                                                                    degree=self._order)

        #gets the total number of control points
        numControlPoints = np.shape(controlPoints_whole)[1] 

        #gets the number of center control points
        numCenterControlPoints =  numControlPoints - 2*self._order

        #gets the number of control points per corridor
        numControlPoints_perCorridor = numIntervalsOfInterestPerCorridor + self._order

        #gets the A and b list
        A_list, b_list = generateConstraintsSFC(sfc_data=sfc_data)


        #calls the function to generate the variables for the controller to create
        controlPoints_cvxpyVar = cp.Variable((2, numControlPoints))

        #now we get the constraints for the control points
        controlPoints_constraints = []

        #iterates over all of the A and b matrices
        for i, (A_temp, b_temp) in enumerate(zip(A_list, b_list)):

            #gets the start index
            startIndex = i * numIntervalsOfInterestPerCorridor
            #gets the end index
            endIndex = startIndex + numControlPoints_perCorridor

            #gets the current section of control points variables
            tempControlPoints_var = controlPoints_cvxpyVar[:,startIndex:endIndex]

            #creates a tiled b temp 
            b_temp_tiled = np.tile(b_temp, (1, numControlPoints_perCorridor))

            #now we set the constraints on the temp control points vector
            tempConstraint = [A_temp @ tempControlPoints_var <= b_temp_tiled]

            #adds the temp constraint to the list of constraints
            controlPoints_constraints += tempConstraint
            
        #now we need to set the equality constraints for the beginning and end variables
        #which will allow us to set them to the values of the start and end control points
        
        #creates the start equality constraint
        startControlPoints_var = controlPoints_cvxpyVar[:,0:self._order]

        #creates the start equality constraint
        startEqualityContraint = [startControlPoints_var == wholePath_startControlPoints]

        #adds this equality constraint to the whole constraints
        controlPoints_constraints += startEqualityContraint

        #creates the end equality constraint
        endControlPoints_var = controlPoints_cvxpyVar[:,-(self._order):]

        #creates the contstraint
        endEqualityConstraint = [endControlPoints_var == wholePath_endControlPoints]

        #adds this back in
        controlPoints_constraints += endEqualityConstraint

        #startEqualityConstraint = [startControlPoints_var = 0]



        #gets the velocity control points
        velocityControlPoints_cp = controlPoints_cvxpyVar[:,0:-1] - controlPoints_cvxpyVar[:,1:]


        #creates the objective function to adjust the control points for this thing.
        minimizeLength_objectiveFunction = cp.Minimize(cp.sum(cp.norm(velocityControlPoints_cp, axis=1)))

        #creates the problem function, which we can then minimize
        prob = cp.Problem(minimizeLength_objectiveFunction, controlPoints_constraints)

        #calls the function to solve the problem 
        prob.solve(solver=cp.CLARABEL)


        print("Status: ", prob.status)
        print("objective: ", prob.value)
        print("Optimized central control points: ", controlPoints_cvxpyVar.value)
        
        #gets the number of points
        potato = 0


    #defines the function to generate the a path by modifying the start and end positions and the theta for velocity
    
    #defines the function to generate control points under a test situation
    def generatePath_test_old(self,
                          numIntervalsOfInterestPerCorridor: int,
                          initialControlPoints: MSG_Control_Points,
                          sfc_data: SFC_Data = None):
        
        #gets the start, center, and end control points
        controlPoints_parsed = initialControlPoints.getControlPointsArray_parsed()

        controlPoints_start = controlPoints_parsed[0]
        controlPoints_center = controlPoints_parsed[1]
        controlPoints_end = controlPoints_parsed[2]

        controlPoints_start = controlPoints_start[:-1,:]
        controlPoints_center = controlPoints_center[:-1,:]
        controlPoints_end = controlPoints_end[:-1,:]


        #gets the main sfc
        sfc = (sfc_data.get_sfc_list())[0]

        #gets the normal vectors
        vertices, normalVectors = sfc.getNormalsVertices_2d()

        normalVectors_old, vertices_old = sfc.getNormalsVertices_old()

        numControlPoints = numIntervalsOfInterestPerCorridor + self._order
        

        #gets the A matrix and b matrix for this
        A, b = generate_A_b(normalVectors=normalVectors_old,
                            vertices=vertices_old)
        
        #gets the tiled version of b
        #TODO figure out whether tiling b is a legitimate thing or not
        b_tiled = np.tile(b, (1, numControlPoints))

        #calls the function to generate the variables for the controller to create
        controlPoints_cvxpyVar = cp.Variable((2, numControlPoints))

        #creates the control points constraints
        controlPoints_constraints = []

        tempConstraint = [A @ controlPoints_cvxpyVar <= b]

        controlPoints_constraints += tempConstraint

        #'''
        #then creates the equality constraints
        startControlPoints_var = controlPoints_cvxpyVar[:,:self._order]

        #creates the start equality constraint
        startEqualityConstraint = [startControlPoints_var == controlPoints_start]

        #adds this into the constraints
        controlPoints_constraints += startEqualityConstraint


        endControlPoints_var = controlPoints_cvxpyVar[:,(-self._order):]

        endEqualityConstraint = [endControlPoints_var == controlPoints_end]

        controlPoints_constraints += endEqualityConstraint
        #'''
        


        velocityControlPoints_cp = controlPoints_cvxpyVar[:,0:-1] - controlPoints_cvxpyVar[:,1:]

        minimizeLength_objectiveFunction = cp.Minimize(cp.sum(cp.norm(velocityControlPoints_cp, axis=1)))

        problem = cp.Problem(minimizeLength_objectiveFunction, controlPoints_constraints)

        problem.solve(solver=cp.CLARABEL)


        print("Status: ", problem.status)
        print("Objective: ", problem.value)
        print("Optimized central control points: ", controlPoints_cvxpyVar.value)


        #gets A times the two control points sets
        startControlPoints_projection = A @ controlPoints_start
        endControlPoints_projection = A @ controlPoints_end

        print("start projection: ", startControlPoints_projection)
        print("end projection: ", endControlPoints_projection)
        print("b: ", b)

        return 0



    #creates the new generate path test function
    def generatePath_test(self,
                          startControlPoints: np.ndarray,
                          endControlPoints: np.ndarray,
                          M: int,
                          sfc_data: SFC_Data = None):
        

        sfc_list = sfc_data.get_sfc_list()

        #gets the start sfc
        sfc_start = sfc_list[0]

        


        potato = 0

    #creates the function to generate the a new path for a singl flight corridor
    def generatePath_test_single(self,
                                 startControlPoints: np.ndarray,
                                 endControlPoints: np.ndarray,
                                 M: int,
                                 sfc: SFC = None):
        


        potato = 0


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
        vertices_temp, normalVectors_temp = sfc_temp.getNormalsVertices_2d()

        #gets the A and b matrices
        A_temp, b_temp = generate_A_b(normalVectors=normalVectors_temp,
                                      vertices=vertices_temp)
        
        A_list.append(A_temp)

        b_list.append(b_temp)

    #returns the A and b lists
    return A_list, b_list



#defines the function to expand the initial and final safe flight corridors
#by incremental steps to include the first and last control points, which have
#been giving me some errors as of late. 'Tis rather frustrating at present, I say.

def expandStartEndSFCs(startSFC: SFC,
                       endSFC: SFC,
                       startControlPoints: np.ndarray,
                       endControlPoints: np.ndarray):
    

    potato = 0