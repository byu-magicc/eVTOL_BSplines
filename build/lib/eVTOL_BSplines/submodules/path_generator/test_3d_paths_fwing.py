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
velocity_scale = 300
start_point = Waypoint(location = np.array([600,-100])[:,None],
                         velocity = np.array([-26,  0])[:,None])
end_point = Waypoint(location=np.array([600,  400])[:,None],
                      velocity=np.array([-26, 0])[:,None])
waypoint_a = Waypoint(location = np.array([600,-100])[:,None], velocity = np.array([-26,  0])[:,None]*velocity_scale)
waypoint_b = Waypoint(location=np.array([600,  400])[:,None], velocity=np.array([-26, 0])[:,None]*velocity_scale)
waypoint_data = WaypointData(start_waypoint=waypoint_a, end_waypoint=waypoint_b)


objective_function_type="minimal_distance_path"
objective_function_type="minimal_velocity_path"
# objective_function_type="minimal_acceleration_path"
max_curvature= 0.0056
# max_curvature=None

fig = plt.figure()

# =============
# First subplot
# =============
# set up the axes for the first plot
order = 3
dimension = 2
num_intervals_free_space = 32
path_gen = PathGenerator(dimension,num_intervals_free_space=num_intervals_free_space)
start_time = time.time()
control_points, status = path_gen.generate_path(waypoint_data=waypoint_data, max_curvature=max_curvature,
    max_incline=None, sfc_data=None, obstacles=None,objective_function_type=objective_function_type,
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


ax = fig.add_subplot()
ax.plot(spline_data[0,:], spline_data[1,:])
ax.set_ylabel("y (m)")
ax.set_xlabel("x (m)")
#  \n \n evaluation time: " + str(np.round(eval_time,2))
# ax.set_aspect('equal')
plot2D_waypoints(waypoint_data, ax, arrow_scale=0.2)

plt.show()

bspline_1 = BsplineEvaluation(control_points, order, 0, 2)
bspline.plot_derivative_magnitude(1000,1)

