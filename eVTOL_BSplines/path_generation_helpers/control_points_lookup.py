#does the same thing as the control points matrix file, but accomplishes this
#using lookup tables instead of recursions
import numpy as np
import os, sys
from pathlib import Path
#imports the reader and the generators
from eVTOL_BSplines.path_generation_helpers.lookUpHelpers import lookUpTableGenerator, lookUpTableReader

#lookup table use for control points generator
class controlPointsGeneratorLookup:

    #creates the initialization function
    def __init__(self,
                 M: int, #the number of intervals of interest for the problem
                 degree: int, #the degree of the polynomial for the problem
                 dimension: int): #the number of dimensions for the problem
        
        #instantiates the lookup generator, which automatically reads up the Lookup tables
        self.reader = lookUpTableReader()


        #saves the variables
        self.M = M
        self.degree = degree
        self.dimension = dimension
      

        tomatoes = 0



    #creates the function to generate the control points:

    def generateControlPoints(self,
                              Start: np.ndarray,
                              End: np.ndarray,
                              rho: np.ndarray):
        

        #values shape checking section:
        #####################################################################################
        #checks to make sure they are the correct sizes
        StartShape = np.shape(Start)
        #raises a ValueError if it is not the right shape for the start
        if StartShape.item(0) != self.dimension:
            raise ValueError(f"Start Elements Dimension: {StartShape[0]} is not equal to set dimension: {self.dimension}")
        EndShape = np.shape(End)
        if StartShape[0] != self.dimension:
            raise ValueError(f"End Elements Dimension: {EndShape[0]} is not equal to set dimension: {self.dimension}")
        #sets the rho length
        rhoLength = np.size(rho)
        if rhoLength != self.degree:
            raise ValueError(f"rho Length: {rhoLength} is not equal to polynomial degree: {self.degree}")

        #####################################################################################

        #concatenates together the boundary conditions
        boundaryConditions = np.concatenate((Start, End), axis=1)

        #gets the solution for the inverse portion
        InversePortion = self.solveInversePortion(rho=rho)


        
    


    #defines the function to obtain the solution for the matrix inverse problem
    def solveInversePortion(self, rho: np.ndarray):
        

        #calls the function to get Y Array sum
        Y_sum = self.sumYArrays(rho=rho)

        #calls the function to obtain the Z array sum
        Z_sum = self.sumZArrays(rho=rho)

        #gets the transposed version of the Y and Z sums
        Y_sum_T = np.transpose(Y_sum)


        #z sum is the one on the same side of the equation as X
        Z_sum_T = np.transpose(Z_sum)

        #gets the X_transpose solution
        X_T = np.linalg.lstsq(Z_sum_T, Y_sum_T, rcond=None)

        #gets the transpose of the X_transpose solution to get the X solution
        X = np.transpose(X_T)

        #returns the X array
        return X


    #defines the function to sum up the y arrays
    def sumYArrays(self, rho: np.ndarray):
        

        #creates the shape 
        Y_sum_height = self.M + self.degree
        Y_sum_width = Y_sum_height - 2*self.degree

        Y_sum = np.zeros((self.M + self.degree, Y_sum_width))
        #iterates over all of the Y arrays corresponding to the rhos
        for ell in range(1, (len(rho)+1)):
            
            
            #creates the Y key
            Y_key = f'degree{self.degree}_l{ell}_M{self.M}'
            Y_temp = (self.reader.Y_data)[Y_key]

            #adds the Y_temp to the Y_sum
            Y_sum = Y_sum + Y_temp

        #returns the Y sum
        return Y_sum
    
    #defines the function to sum up the z arrays
    def sumZArrays(self, rho: np.ndarray):

        Z_sum_width = self.M - self.degree
        Z_sum_height = Z_sum_width

        Z_sum = np.zeros((Z_sum_height, Z_sum_width))

        for ell in range(1, (len(rho) + 1)):

            #creates the Y key
            Z_key = f'degree{self.degree}_l{ell}_M{self.M}'
            Z_temp = (self.reader.Z_data)[Z_key]

            #adds the Z_temp to the Z_sum
            Z_sum = Z_sum + Z_temp
        
        #returns the Z sum
        return Z_sum
    
    #defines the function to 