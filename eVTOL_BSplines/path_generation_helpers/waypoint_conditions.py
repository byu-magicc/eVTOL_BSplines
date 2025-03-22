#this file implements the class for the waypoint conditions
import numpy as np


#creates the waypoint conditions class
class waypoint_conditions:

    #creates the initialization function
    def __init__(self,
                 position: np.ndarray,
                 velocity: np.ndarray,
                 acceleration: np.ndarray,
                 jerk: np.ndarray,
                 snap: np.ndarray):
        
        #saves the position velocity and accel vectors
        self.position = position
        self.velocity = velocity
        self.acceleration = acceleration
        self.jerk = jerk
        self.snap = snap
