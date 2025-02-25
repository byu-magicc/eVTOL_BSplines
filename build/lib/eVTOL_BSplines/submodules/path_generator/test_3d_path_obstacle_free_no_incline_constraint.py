import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from bsplinegenerator.bsplines import BsplineEvaluation
# from path_generation.path_generator import PathGenerator
from path_generation.path_generator import PathGenerator
from path_generation.safe_flight_corridor import SFC, SFC_Data, plot_sfcs, get2DRotationAndTranslationFromPoints
from path_generation.obstacle import Obstacle, plot_2D_obstacles, plot_3D_obstacles, plot_cylinders
from path_generation.waypoint_data import plot2D_waypoints, plot3D_waypoints
from path_generation.waypoint_data import Waypoint, WaypointData
from path_generation.path_plotter import set_axes_equal
import time


# waypoints
velocity_scale = 100
waypoint_a = Waypoint(location=np.array([[0],[0],[-300]]), velocity = np.array([[0.85],[0.5267826876426369],[0]])*velocity_scale)
# waypoint_b = Waypoint(location=np.array([[-100],[-50],[0]]), velocity = np.array([[1],[0],[0]])*velocity_scale)
waypoint_b = Waypoint(location=np.array([[600],[50],[0]]), velocity = np.array([[1],[0],[0]])*velocity_scale)
waypoint_data = WaypointData(start_waypoint=waypoint_a, end_waypoint=waypoint_b)



#obstacle_list

obstacles = None

objective_function_type="minimal_distance_path"
max_curvature= 0.0114362
# max_incline = 0.2679491924311227
max_incline = None
max_curvature=None

fig = plt.figure()

# =============
# First subplot
# =============
# set up the axes for the first plot
order = 3
dimension = 3
num_intervals_free_space = 5
path_gen = PathGenerator(dimension, num_intervals_free_space = num_intervals_free_space)
start_time = time.time()
control_points, status = path_gen.generate_path(waypoint_data=waypoint_data, max_curvature=max_curvature,
    max_incline=max_incline, sfc_data=None, obstacles=obstacles,objective_function_type=objective_function_type,
    obstacle_type="cylinder")
end_time = time.time()
eval_time = end_time - start_time
num_cont_pts = np.shape(control_points)[1]
print("control_points: " , control_points)
print("computation time:" , end_time - start_time)

spline_start_time = 0
scale_factor = 1
bspline = BsplineEvaluation(control_points, order, spline_start_time, scale_factor, False)
curvature_data, time_data = bspline.get_spline_curvature_data(10000)
velocity_data, time_data = bspline.get_spline_derivative_data(10000, 1)
z_vel_data = velocity_data[2,:]
horiz_data = np.linalg.norm(velocity_data[0:2,:],2,0)
incline_data = z_vel_data / horiz_data

number_data_points = 10000
spline_data, time_data = bspline.get_spline_data(number_data_points)

path_length = bspline.get_arc_length(10000)

# distances_to_obst = np.linalg.norm(spline_data[0:2,:] - obstacle_center[0:2,:], 2, 0)
# closest_distance_to_obst = np.min(distances_to_obst) - obstacle_radius
# print("dist to obst: " , closest_distance_to_obst)
print("path length: " , path_length)
print("max curvature: " , np.max(curvature_data))
print("max incline: " , np.max(np.abs(incline_data)))

ax = fig.add_subplot(projection='3d')
ax.plot(spline_data[0,:], spline_data[1,:], spline_data[2,:])
ax.set_ylabel("y")
ax.set_xlabel("x")
ax.set_zlabel("z")
#  \n \n evaluation time: " + str(np.round(eval_time,2))
# ax.set_aspect('equal')
plot3D_waypoints(waypoint_data, ax, arrow_scale=3)
# plot_cylinders(obstacles, ax)
set_axes_equal(ax,dimension)
# plt.axis('off')
plt.show()

