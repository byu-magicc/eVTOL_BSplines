import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from bsplinegenerator.bsplines import BsplineEvaluation
from path_generation.path_generator import PathGenerator
from path_generation.safe_flight_corridor import SFC_Data, get3DRotationAndTranslationFromPoints
from path_generation.path_plotter import set_axes_equal
from path_generation.waypoint_data import Waypoint, WaypointData
import time


waypoint_1 = Waypoint(location=np.array([[0],[0],[100]]))
waypoint_2 = Waypoint(location=np.array([[0],[100],[105]]))
waypoint_1.velocity = np.array([[1],[0],[0]])
waypoint_2.velocity = np.array([[-1],[ 0],[ 0]])
waypoint_data = WaypointData(start_waypoint=waypoint_1,end_waypoint=waypoint_2)

dimension = np.shape(waypoint_1.location)[0]
max_curvature = 1
# max_curvature = None
max_incline = None
order = 3


path_gen = PathGenerator(dimension)
start_time = time.time()
# objective_function_type = "minimal_distance_path"
# objective_function_type = "minimal_velocity_path"
objective_function_type = "minimal_acceleration_path"
control_points, status = path_gen.generate_path(waypoint_data=waypoint_data, max_curvature=max_curvature,
    max_incline=max_incline, sfc_data=None, obstacles=None, objective_function_type=objective_function_type)
end_time = time.time()
print("control_points: " , control_points)

# print("control_points: " , control_points)
spline_start_time = 0
scale_factor = 1
bspline = BsplineEvaluation(control_points, order, spline_start_time, scale_factor, False)
number_data_points = 10000
spline_data, time_data = bspline.get_spline_data(number_data_points)
spline_velocity_data, time_data = bspline.get_spline_derivative_data(number_data_points,1)
z_vel_mag = spline_velocity_data[2,:]
horiz_vel_mag = np.linalg.norm(spline_velocity_data[0:2,:],2,0)
incline_data = z_vel_mag/horiz_vel_mag
curvature_data, time_data = bspline.get_spline_curvature_data(number_data_points)
velocity_data, time_data = bspline.get_derivative_magnitude_data(number_data_points,1)
# acceleration_data, time_data = bspline.get_spline_derivative_data(1000,2)
minvo_cps = bspline.get_minvo_control_points()
waypoints = waypoint_data.get_waypoint_locations()
end_time_spline = bspline.get_end_time()
start_velocity = bspline.get_derivative_at_time_t(0,1)
start_acceleration = bspline.get_derivative_at_time_t(0,2)
end_velocity = bspline.get_derivative_at_time_t(end_time_spline,1)
end_acceleration = bspline.get_derivative_at_time_t(end_time_spline,2)
path_length = bspline.get_arc_length(number_data_points)
print("start_velocity: " , start_velocity)
print("start_acceleration: " , start_acceleration)
print("end_velocity: " , end_velocity)
print("end_acceleration: " , end_acceleration)
print("max incline" , np.max(incline_data))
print("max curvature" , np.max(curvature_data))
print("start curvature: " , curvature_data[0])
print("end curvature: " , curvature_data[-1])
print("path length: " , path_length)
print("computation time: " , end_time - start_time)
# print("end curvature: " , curvature_data[-1])

# print("sfcs: " , sfcs)
plt.figure()
ax = plt.axes(projection='3d')
ax.plot(spline_data[0,:], spline_data[1,:],spline_data[2,:])
ax.scatter(control_points[0,:], control_points[1,:],control_points[2,:])
# ax.scatter(waypoints[0,:],waypoints[1,:],waypoints[2,:])
# plt.scatter(minvo_cps[0,:],minvo_cps[1,:])
ax.set_xlabel('$X$')
ax.set_ylabel('$Y$')
ax.set_zlabel('$Z$')
set_axes_equal(ax,dimension)
plt.show()

# plt.figure()
# plt.plot(time_data, incline_data)
# plt.show()

# plt.figure()
# plt.plot(time_data, curvature_data)
# plt.show()

# plt.figure()
# plt.plot(time_data, velocity_data)
# plt.show()


