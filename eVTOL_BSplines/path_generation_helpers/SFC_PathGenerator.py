#this is the file which implements the ability to generate 
#bspline paths through safe flight corridors. But this will
#be implemented differently from David Christensen's work, because
#it will be done by 

import numpy as np

class SFC_PathGenerator:

    #creates the init function
    def __init__(self,
                 dimension: int = 2,
                 num_intervals_free_space: int = 5):
        
        #sets the degree being used here
        self.degree = 3

        #saves the arguments
        self.dimension = dimension
        self.num_intervals_free_space = num_intervals_free_space
        pass

    #creates the function to generate the path
    def generate_path(self,
                      start_conditions: list[np.ndarray],
                      end_conditions: list[np.ndarray]):
