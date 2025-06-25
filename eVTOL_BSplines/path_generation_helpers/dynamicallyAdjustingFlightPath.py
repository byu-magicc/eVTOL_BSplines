#This file implements a dynamically adjusting flight path planning class
#for the eVTOL problem

import numpy as np

#instantiates the lookupTable reader class
from eVTOL_BSplines.path_generation_helpers.lookup_table_helpers import lookUpTableReader, YZGeneratorReader


#creates the dynamic flight path, which does not stay constant, but increments and updates 
class dynamicFlightPath:

    #creates the initialization function
    #Arguments:
    #1. the initial conditions for the whole spline, at the actual start
    #2. the final conditions for the whole spline
    #3. the number of dimensions for this problem
    #4. the degree of the bspline for this problem
    #5. the number of intervals of interest for this
    def __init__(self,
                 initialConditionsMain: list[np.ndarray],
                 finalConditionsMain: list[np.ndarray],
                 numDimensions: int = 2,
                 d: int = 3,
                 M: int = 100):
        
        #saves all of those arguments here
        self.initialConditionsMain = initialConditionsMain
        self.finalConditionsMain = finalConditionsMain
        self.numDimensions = numDimensions
        self.d = d
        self.M = M


        #creates the reader which loads it automatically
        self.reader = lookUpTableReader()

        #creates the YZ helper and loads it all up
        self.YZ_helper = YZGeneratorReader()
        (self.YZ_helper).loadYZTables()


        #creates the Current M variable, which starts off as the initial M
        #but of course it will change and shrink until M = d
        self.current_M = self.M
        pass


    #defines the function that iterates over all of the control points
    #def updatePath(self):






    #defines the function to get the current set of control points
    def getCurrentControlPoints(self,
                                current_M: int, #puts in the current number of initial conditions, which gets progressively smaller
                                rho: np.ndarray, #the mixing matrix to mix importance of each derivative
                                currentInitialConditions: list[np.ndarray], #the current conditions of the aircraft
                                currentFinalConditions: list[np.ndarray]): #the final conditions of the aircraft (Which I think should remain constant and not update. But I don't know what I'm going to want)

        #concatenates together the current start and end conditions
        completeCurrentConditions = self.combineStartEndConditions(startConditions=currentInitialConditions,
                                                                   endConditions=currentFinalConditions)
        
        #gets the pseudoinverse we need
        B_cat_pseudoinverse = self.reader.getIndividualPseudoinverse(M=current_M,
                                                                     d=self.d)
        
        #gets the U2 we need
        U2 = self.reader.getIndividualU2(M=current_M,
                                         d=self.d)

        #gets the inverted portion
        X = self.getInvertedPortion(current_M=current_M,
                                    rho=rho)

        #gets the control points
        CtrlPoints = completeCurrentConditions @ B_cat_pseudoinverse @ (np.eye(current_M + self.d) - X @ np.transpose(U2))

        #and returns them back
        return CtrlPoints


    #creates the function to get a particular number of control points for the current conditions
    def getInvertedPortion(self,
                           current_M: int,
                           rho: np.ndarray):
        
        #creates the matrices to store the  summation matrices 
        Y_sum = np.zeros(((current_M + self.d), (current_M - self.d)))
        Z_sum = np.zeros(((current_M - self.d, current_M - self.d)))

        #creates the counter variable l, which starts at one and goes up to degree inclusive
        l = 1

        #iterates over all of the items in rho, to ge the summations for Z and Y
        for rho_l in rho:

            #gets the temporary Y and the temporary Z variables
            temp_Y = self.YZ_helper.getIndividualYMatrix(M=current_M,
                                                         d=self.d,
                                                         l=l)
            
            #gets the  same for the Z
            temp_Z = self.YZ_helper.getIndividualZMatrix(M=current_M,
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