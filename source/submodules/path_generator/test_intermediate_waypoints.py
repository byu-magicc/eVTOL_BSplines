import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from bsplinegenerator.bsplines import BsplineEvaluation
# from path_generation.path_generator import PathGenerator
from path_generation.path_generator import PathGenerator
from path_generation.safe_flight_corridor import SFC_Data, get3DRotationAndTranslationFromPoints
from path_generation.path_plotter import set_axes_equal
from path_generation.waypoint_data import Waypoint, WaypointData, plot3D_waypoints
from path_generation.obstacle import Obstacle, plot_3D_obstacles
import time

#note incline constraints work much better when have a start and an end direction

dimension = 3
max_curvature = 1
max_incline = 2
# max_incline = None
order = 3
scale_factor = 1
# path_objective_type = "minimal_acceleration_path"
path_objective_type = "minimal_velocity_path"

##### Initial Settings
# path 1
waypoint_1_loc = np.array([[3],[4],[0]])
waypoint_2_loc = np.array([[7],[10],[13]])
waypoint_2_5_loc = np.array([[ 9.80919643],[ 1.8519873 ],[10.37864785]])
waypoint_3_loc = np.array([[8],[-10],[7]])
waypoint_sequence = np.concatenate((waypoint_1_loc, waypoint_2_loc, waypoint_2_5_loc, waypoint_3_loc),1)
start_waypoint = Waypoint(location=waypoint_1_loc)
end_waypoint= Waypoint(location=waypoint_3_loc)
waypoint_data = WaypointData(start_waypoint, end_waypoint)
waypoint_data.set_location_data(waypoint_sequence)

obstacles = [Obstacle(np.array([[4.5],[5.5],[4]]), 2.0),
               Obstacle(np.array([[10],[-5],[8]]), 1.5)]

path_gen = PathGenerator(dimension, num_intervals_free_space=8)
start_time = time.time()
control_points, status = path_gen.generate_path(waypoint_data=waypoint_data, max_curvature=max_curvature,
    max_incline=max_incline, sfc_data=None, obstacles=obstacles)
end_time = time.time()
print("computation time: " , end_time - start_time)

spline_start_time = 0
bspline = BsplineEvaluation(control_points, order, spline_start_time, scale_factor, False)
number_data_points = 10000
spline_data, time_data = bspline.get_spline_data(number_data_points)
curvature_data, time_data  = bspline.get_spline_curvature_data(number_data_points)


plt.figure()
ax = plt.axes(projection='3d')
ax.plot(spline_data[0,:], spline_data[1,:],spline_data[2,:])
plot3D_waypoints(waypoint_data,ax)
plot_3D_obstacles(obstacles,ax)
ax.set_xlabel('$X$')
ax.set_ylabel('$Y$')
ax.set_zlabel('$Z$')
set_axes_equal(ax,dimension)
plt.show()

plt.figure()
plt.title("Curvature")
plt.plot(time_data, curvature_data,color = "b")
plt.show()
