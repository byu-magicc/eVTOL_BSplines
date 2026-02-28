#This file implements the constant or static flight path for the conditions
import numpy as np

#instantiates the lookupTable reader class
from eVTOL_BSplines.path_generation_helpers.lookup_table_helpers import lookUpTableReader, YZGeneratorReader

#creates the static flight path, which stays constant from start to finish, and does not update recursively
class staticFlightPath:

    #saves the main initial conditions, but which will be updated as we move along
    def __init__(self):

        #instantiates the reader
        self.reader = lookUpTableReader()
        self.YZ_helper = YZGeneratorReader()

        #calls the function to load the Y and Z tables
        (self.YZ_helper).loadYZTables()
        pass


    #solves for the Control points
    def getControlPoints(self,
                         initialConditions: list[np.ndarray],
                         finalConditions: list[np.ndarray],
                         rho: np.ndarray,
                         numDimensions: int = 2,
                         d: int = 3,
                         M: int = 10):

        #saves the numbers
        self.numDimensions = numDimensions
        self.d = d
        self.M = M

        #gets the conditions matrix
        completeConditions = self.combineStartEndConditions(startConditions=initialConditions,
                                                            endConditions=finalConditions)
        
        #gets the pseudoinverse portion that we need
        B_cat_pseudoinverse = self.reader.getIndividualPseudoinverse(M=self.M, d=self.d)

        #gets the U2
        U2 = self.reader.getIndividualU2(M=self.M, 
                                         d=self.d)
        
        #gets the inverted portion
        X = self.getInvertedPortion(rho=rho)
        
        #gets the control points
        CtrlPoints = completeConditions @ B_cat_pseudoinverse @ (np.eye(self.M + self.d) - X @ np.transpose(U2))

        #returns the control points
        return CtrlPoints
    
    #defines the function to get the localized control points
    def getLocalizedControlPoints(self,
                                  conditions: list[np.ndarray],
                                  d: int,
                                  M: int):

        #calls the reader function to get the individual B inverse matrix
        B_hat_inv = self.reader.getIndividualBHatInv(M=M, d=d)

        #concatenates together the conditions
        conditionsArray = np.concatenate(conditions, axis=1)

        ctrlPts_localized = conditionsArray @ B_hat_inv

        #returns the control points
        return ctrlPts_localized

    #from d control points, get the corresponding condition
    #(Go the opposite way as the previous function)
    def getLocalizedConditions(self,
                               controlPoints: np.ndarray,
                               d: int,
                               M: int):

        #calls the reader function to get the individual B inverse matrix
        B_hat_inv = self.reader.getIndividualBHatInv(M=M, d=d)

        #gets the inverse of Bhat inverse. Because I'm lazy
        B_hat = np.linalg.inv(B_hat_inv)

        #gets the startConditions
        localizedConditions = controlPoints @ B_hat

        return localizedConditions



    
    #creates the function to get the inverted portion, which helps break up the problem a little bit
    def getInvertedPortion(self,
                           rho: np.ndarray): #takes in the rho scaling vector for this portion
        
        #creates the matrices to store the  summation matrices 
        Y_sum = np.zeros(((self.M + self.d), (self.M - self.d)))
        Z_sum = np.zeros(((self.M - self.d, self.M - self.d)))

        #creates the counter variable l, which starts at one and goes up to degree inclusive
        l = 1

        #iterates over all of the items in rho, to ge the summations for Z and Y
        for rho_l in rho:

            #gets the temporary Y and the temporary Z variables
            temp_Y = self.YZ_helper.getIndividualYMatrix(M=self.M,
                                                         d=self.d,
                                                         l=l)
            
            #gets the  same for the Z
            temp_Z = self.YZ_helper.getIndividualZMatrix(M=self.M,
                                                         d=self.d,
                                                         l=l)
            
            #adds it to the rho scaled sum
            Y_sum = Y_sum + rho_l*temp_Y

            #same for the Z sum
            Z_sum = Z_sum + rho_l*temp_Z


            #increments l by one
            l = l + 1


        #after obtaining the sum, we need to explain this. This is important to explain what is happening for posterity

        #so in the main formula for the Control points, it looks something like:
        #C* = [S E] * B_cat_pseudoinverse * (I - W(rho) U2 (U2^T W(rho) U2)^-1 U2^T))
        #And I need an efficient way to calculate that inner (U2^T W(rho) U2)^-1.
        #but I rememeber that one should avoid taking the direct inverse as much as humanly possible.

        #defines the Z_sum bar
        Z_sum_bar = np.transpose(Z_sum)

        #defines the Y_sum_bar
        Y_sum_bar = np.transpose(Y_sum)

        #solves for the X variable using the linear algebra solve function
        X_bar = np.linalg.solve(Z_sum_bar, Y_sum_bar)

        #gets X
        X = np.transpose(X_bar)

        #returns the X
        return X
    

    #creates the function that get the conditions as a matrix, from a list
    def conditionsMatrixFromList(self, 
                                 conditions: list[np.ndarray]):
        
        #creates the conditions matrix, on which to concatenate everything.
        conditionsMatrix = np.ndarray((self.numDimensions,0))
        #iterates through all the conditions
        for temp_condition in conditions:
            
            conditionsMatrix = np.concatenate((conditionsMatrix, temp_condition), axis=1)


        #returns the whole put together conditions matrix
        return conditionsMatrix
    
    #creates a function to concatenate together the start and end conditions
    def combineStartEndConditions(self,
                                  startConditions: list[np.ndarray],
                                  endConditions: list[np.ndarray]):
        
        #gets the start and end conditiosn matrices
        startMatrix = self.conditionsMatrixFromList(startConditions)
        endMatrix = self.conditionsMatrixFromList(endConditions)

        #concatentates thrm together
        completeConditions = np.concatenate((startMatrix, endMatrix), axis=1)

        return completeConditions
    

    #function to get control points from a list of sequential conditions. That is you might have a list
    #of 5 points each with their own sets of position, velocity, and acceleration conditions. And I want to 
    #get the overlapping sequences of control points from that.
    def getControlPoints_conditionsList(self,
                                        conditions_list: list[list[np.ndarray]],
                                        rho: np.ndarray,
                                        numDimensions: int = 2,
                                        d: int = 3,
                                        M: int = 10):
        controlPoints = np.ndarray((numDimensions,0))
        #iterates over all the conditions in the list
        numConditions = len(conditions_list)
        for i in range(numConditions - 1):

            startConditions = conditions_list[i]
            endConditions = conditions_list[i+1]
            #gets the control points using the big function
            tempControlPoints = self.getControlPoints(initialConditions=startConditions,
                                                      finalConditions=endConditions,
                                                      rho=rho,
                                                      numDimensions=numDimensions,
                                                      d=d,
                                                      M=M)
            
            if i == 0:
                #gets the partition of the control points
                tempControlPoints_section = tempControlPoints
            else:
                tempControlPoints_section = tempControlPoints[:,d:]

            controlPoints = np.concatenate((controlPoints, tempControlPoints_section), axis=1)

        return controlPoints
            
