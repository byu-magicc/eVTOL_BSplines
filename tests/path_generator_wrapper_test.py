import os, sys
from pathlib import Path
import numpy as np
import time
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from bsplinegenerator.bsplines import BsplineEvaluation
tempPath = sys.path
#inserts the parent directory
sys.path.insert(0,os.fspath(Path(__file__).parents[1]))
tempPath = sys.path

from wrappers.path_generator_simplified import path_generator_simplified

from wrappers.path_generator_simplified import Waypoint, WaypointData

pathGen = path_generator_simplified(dimension=2)

#creates the waypoints
startWaypoint = Waypoint(location=np.array([[0.0],[0.0]]), velocity=np.array([[1.0],[0.0]]))
endWaypoint = Waypoint(location=np.array([[10.0],[5.0]]), velocity=np.array([[0.0],[-2.0]]))

waypoints = WaypointData(start_waypoint=startWaypoint,end_waypoint=endWaypoint)

#creates the instance of the path generator
pathGen = path_generator_simplified(dimension=2)

#gets the data
splineData, contrlPoints, timeData = pathGen.waypointPathGenerator(waypoints=waypoints,
                                                                   max_curvature=100.0,
                                                                   degree=3,
                                                                   num_points_per_interval=100)



