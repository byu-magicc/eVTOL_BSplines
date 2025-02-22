import os, sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from bsplinegenerator.bsplines import BsplineEvaluation
# from path_generation.path_generator import PathGenerator
from path_generation.path_generator import PathGenerator
from path_generation.safe_flight_corridor import SFC, SFC_Data, plot_sfcs, get2DRotationAndTranslationFromPoints
from path_generation.obstacle import Obstacle, plot_2D_obstacles, plot_3D_obstacles, plot_cylinders, plot_2d_buildings, \
 plot_buildings
from path_generation.waypoint_data import plot2D_waypoints, plot3D_waypoints
from path_generation.waypoint_data import Waypoint, WaypointData
from path_generation.path_plotter import set_axes_equal
import time

obstaclesPath = os.path.abspath(os.path.join(os.path.dirname(__file__), 'obstacles.npy'))
obstacle_data = np.load(obstaclesPath)

obstacle = Obstacle(np.array([250,250,0]), 282.84/2, 300)
obstacle_list = []
num_obstacles = np.shape(obstacle_data)[1]
for i in range(num_obstacles):
    center_ = np.array([obstacle_data[0,i], obstacle_data[1,i], 0])
    height_ = obstacle_data[2,i]
    radius_ = obstacle_data[3,i]
    obstacle_ = Obstacle(center= center_, radius=radius_, height=height_)
    obstacle_list.append(obstacle_)

# obstacles = []
sfc1 = SFC(dimensions = np.array([[133.33], [633.33], [300]]),
           translation = np.array([[0], [250], [-150]]),
           rotation = np.identity(3))
sfc2 = SFC(dimensions = np.array([[300], [466.666], [201.2]]),
           translation = np.array([[83.33333], [333.33], [-199.4]]),
           rotation = np.identity(3))
sfc3 = SFC(dimensions = np.array([[633.33], [300], [219.1]]),
           translation = np.array([[250], [416.66666667], [-190.95]]),
           rotation = np.identity(3))

sfc_list = [sfc1, sfc2, sfc3]
# sfc1 = SFC(dimensions = np.array([[133.33], [633.33], [300]]),
#            translation = np.array([[0], [250], [-150]]),
#            rotation = np.identity(3))
# sfc2 = SFC(dimensions = np.array([[633.33], [133.33], [300]]),
#            translation = np.array([[250], [166.66], [-150]]),
#            rotation = np.identity(3))
# sfc3 = SFC(dimensions = np.array([[133.33], [633.33], [300]]),
#            translation = np.array([[166.66], [250], [-150]]),
#            rotation = np.identity(3))
# sfc4 = SFC(dimensions = np.array([[633.33], [133.33], [300]]),
#            translation = np.array([[250], [333.33], [-150]]),
#            rotation = np.identity(3))
# sfc5 = SFC(dimensions = np.array([[133.33], [633.33], [300]]),
#            translation = np.array([[333.33], [250], [-150]]),
#            rotation = np.identity(3))
# sfc6 = SFC(dimensions = np.array([[633.33], [133.33], [300]]),
#            translation = np.array([[250], [500], [-150]]),
#            rotation = np.identity(3))
# sfc_list = [sfc1, sfc2, sfc3, sfc4, sfc5, sfc6]
# sfc_a = SFC(dimensions = np.array([[600], [600], [300]]),
#            translation = np.array([[250], [250], [150]]),
#            rotation = np.identity(3))
# sfc_b = SFC(dimensions = np.array([[300], [300], [300]]),
#            translation = np.array([[400], [400], [150]]),
#            rotation = np.identity(3))
# sfc_list = [sfc_a, sfc_b]

point_sequence = np.array([[0,  0,   83.33333   , 500],
                           [0, 333.33,333.33, 500],
                           [-100,-199.4, -199.4, -200 ]])

sfc_data = SFC_Data(sfc_list, point_sequence, 1)

# waypoints
velocity_scale = 300
waypoint_a = Waypoint(location=np.array([[0],[0],[-100]]), velocity = np.array([[0.85],[0.5267826876426369],[0]])*velocity_scale)
waypoint_b = Waypoint(location=np.array([[500],[500],[-200]]), velocity = np.array([[0],[1],[0]])*velocity_scale)
waypoint_data = WaypointData(start_waypoint=waypoint_a, end_waypoint=waypoint_b)

#obstacle_list

objective_function_type="minimal_distance_path"
max_curvature= 0.0114362
max_incline = 0.17279
# max_curvature=None

fig = plt.figure()

# =============
# First subplot
# =============
# set up the axes for the first plot
order = 3
dimension = 3
num_intervals_free_space = 7 
path_gen = PathGenerator(dimension,num_intervals_free_space=num_intervals_free_space)
start_time = time.time()
control_points, status = path_gen.generate_path(waypoint_data=waypoint_data, max_curvature=max_curvature,
    max_incline=max_incline, sfc_data=sfc_data,objective_function_type=objective_function_type,
    obstacle_type="cylinder")
print("control points: " , control_points)
end_time = time.time()
eval_time = end_time - start_time
num_cont_pts = np.shape(control_points)[1]

spline_start_time = 0
scale_factor = 1
bspline = BsplineEvaluation(control_points, order, spline_start_time, scale_factor, False)
curvature_data, time_data = bspline.get_spline_curvature_data(10000)
path_length = bspline.get_arc_length(100000)
print("max curvature: " , np.max(curvature_data))
print("path length: " , path_length)
print("computation time:" , end_time - start_time)
number_data_points = 10000
spline_data, time_data = bspline.get_spline_data(number_data_points)

ax = fig.add_subplot(projection="3d")
ax.plot(spline_data[0,:], spline_data[1,:], spline_data[2,:] )
# ax.set_ylabel("y (m)")
# ax.set_xlabel("x (m)")
#  \n \n evaluation time: " + str(np.round(eval_time,2))
# ax.set_aspect('equal')
plot3D_waypoints(waypoint_data, ax, arrow_scale=1)
plot_buildings(obstacle_list, ax)
plot_sfcs(sfc_list, ax, 0.5)
plt.axis('off')
plt.show()

