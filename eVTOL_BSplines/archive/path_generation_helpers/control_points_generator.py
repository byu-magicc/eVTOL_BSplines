#this file implements the generation for control points
#for the thing.



import os, sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from bsplinegenerator.bsplines import BsplineEvaluation
from IPython.display import display
import sympy as sp

from IPython.display import display

import scipy.integrate as integrate

sys.path.insert(0,os.fspath(Path(__file__).parents[1]))
tempPath = sys.path



from eVTOL_BSplines.path_generation_helpers.matrix_generators_efficient import create_W_Matrix
from eVTOL_BSplines.path_generation_helpers.matrix_helpers import B_hat_B_hat_inv_simplified


class PathGenerator:


    #creates the initialization function

    def __init__(self,
                 dimension: int, #the dimensionality of the problem
                 M: int, #the number of intervals of interes
                 degree: int = 3, #the degree of the polynomial in question
                 num_intervals_free_space: int=5):

        #saves the dimension
        self.dimension = dimension
        self.M = M
        #saves the degree of the polynomial in question
        self.degree = degree
        self.num_intervals_free_space = num_intervals_free_space

        #instantiates the class to perform W matrix creation
        self.W_matrix_creator = create_W_Matrix(d=degree,integratorFileName = "lookUpTables/degree_5_integrations.npz")


    #creates the function to get the control points
    def getControlPoints(self,
                         S: np.ndarray, #the start conditions for the spline,
                         E: np.ndarray, #the end conditions for the spline,
                         rho: np.ndarray): #the scaling array for the derivatives conditions
        

        #gets the useful submatrices partition
        W_cc, W_tc, W_cr = (self.W_matrix_creator).W_applicable_subsections(d=self.degree,
                                                                            M=self.M,
                                                                            rho=rho)

        #gets the B hat matrix 
        B_hat, B_hat_inv = B_hat_B_hat_inv_simplified(degree=self.degree,
                                                      alpha=1.0)
        
        #gets the first section of control points
        C_S = S @ B_hat_inv

        #gets the final section of control points
        C_E = E @ B_hat_inv

        #now, for the center section, we would like to obtain the control points.
        #To do this, we start with the right side of the equation, which is straightforward to compute
        A_right = -B_hat_inv @ (S @ W_tc.T + E @ W_cr.T)

        

        potato = 0