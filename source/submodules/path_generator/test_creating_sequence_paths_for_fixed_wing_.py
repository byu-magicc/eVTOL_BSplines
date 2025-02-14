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
max_curvature = .009
# max_incline = None
order = 3
scale_factor = 1
# path_objective_type = "minimal_acceleration_path"
path_objective_type = "minimal_velocity_path"

##### Initial Settings
# path 1
waypoint_1 = Waypoint(location=np.array([[0],[0],[-20]]))
waypoint_2 = Waypoint(location=np.array([[0],[350],[-18]]))
waypoint_1.velocity = np.array([[100],[0],[0]])
waypoint_1.acceleration = np.array([[0],[0],[0]])
# waypoint_2.acceleration = np.array([[0],[0],[0]])
waypoint_2.velocity = np.array([[-100],[0],[0]])
waypoint_data_1 = WaypointData(start_waypoint=waypoint_1,end_waypoint=waypoint_2)

# plt.figure()
# ax = plt.axes(projection='3d')
# ax.set_xlabel('$X$')
# ax.set_ylabel('$Y$')
# ax.set_zlabel('$Z$')
# plot3D_waypoints(waypoint_data_1, ax)
# set_axes_equal(ax,dimension)
# plt.show()
number_data_points = 10000

### 1st path
path_gen = PathGenerator(dimension, num_intervals_free_space=5)
start_time_1 = time.time()
control_points, status = path_gen.generate_path(waypoint_data=waypoint_data_1, max_curvature=max_curvature,
    max_incline=None, sfc_data=None, obstacles=None, objective_function_type=path_objective_type)
end_time_1 = time.time()
print("control_points: " , control_points)
print(" ")
print("PATH 1 DATA:")
print("computation time: " , end_time_1 - start_time_1)
spline_start_time_1 = 0
bspline = BsplineEvaluation(control_points, order, spline_start_time_1, scale_factor, False)
end_time_spline_1 = bspline.get_end_time()
#### path 1 data
spline_data_1, time_data_1 = bspline.get_spline_data(number_data_points)
curvature_data_1, time_data_1 = bspline.get_spline_curvature_data(number_data_points)
velocity_data_1, time_data_1 = bspline.get_spline_derivative_data(number_data_points,1)
acceleration_data_1, time_data_2 = bspline.get_spline_derivative_data(number_data_points,2)
z_vel_mag_1 = velocity_data_1[2,:]
horiz_vel_mag_1 = np.linalg.norm(velocity_data_1[0:2,:],2,0)
incline_data_1 = z_vel_mag_1/horiz_vel_mag_1

print("max incline: " , np.max(abs(incline_data_1)))
print("max curvature" , np.max(curvature_data_1))
# print(" ")
# ### plot 1st path
# plt.figure()
# ax = plt.axes(projection='3d')
# ax.plot(spline_data_1[0,:], spline_data_1[1,:],spline_data_1[2,:], color = "b")
# ax.set_xlabel('$X$')
# ax.set_ylabel('$Y$')
# ax.set_zlabel('$Z$')
# plot3D_waypoints(waypoint_data_1, ax)
# set_axes_equal(ax,dimension)
# plt.show()

### second path
waypoint_1_two = Waypoint(location=np.array([[0],[350],[-18]]))
waypoint_2_two = Waypoint(location=np.array([[0],[700],[-22]]))
end_time = bspline.get_end_time()
# waypoint_1_two.acceleration = np.array([0,0,0])
waypoint_1_two.velocity = bspline.get_derivative_at_time_t(end_time,1)
waypoint_1_two.acceleration = bspline.get_derivative_at_time_t(end_time,2)
# print("waypoint_1_two.acceleration: ", waypoint_1_two.acceleration)
waypoint_2_two.velocity = np.array([[100],[0],[0]])
waypoint_data_2 = WaypointData(start_waypoint=waypoint_1_two,end_waypoint=waypoint_2_two)
spline_start_time_2 = end_time_spline_1
start_time_2 = time.time()
# path_gen.set_num_intervals_free_space(7)
control_points_2, status_2 = path_gen.generate_path(waypoint_data=waypoint_data_2, max_curvature=max_curvature,
    max_incline=None, sfc_data=None, obstacles=None, objective_function_type=path_objective_type)
end_time_2 = time.time()
print("control_points_2: " , control_points_2)
## spline 2 data
bspline_2 = BsplineEvaluation(control_points_2, order, spline_start_time_2, scale_factor, False)
spline_data_2, time_data_2 = bspline_2.get_spline_data(number_data_points)
curvature_data_2, time_data_2 = bspline_2.get_spline_curvature_data(number_data_points)
velocity_data_2, time_data_2 = bspline_2.get_spline_derivative_data(number_data_points,1)
acceleration_data_2, time_data_2 = bspline_2.get_spline_derivative_data(number_data_points,2)
z_vel_mag_2 = velocity_data_2[2,:]
horiz_vel_mag_2 = np.linalg.norm(velocity_data_2[0:2,:],2,0)
# print(" z_vel_mag_2 ", z_vel_mag_2)
# print(" horiz_vel_mag_2 ", np.max(horiz_vel_mag_2))
incline_data_2 = z_vel_mag_2/horiz_vel_mag_2

waypoints = np.concatenate((waypoint_data_2.start_waypoint.location, waypoint_data_2.end_waypoint.location),1)
print("waypoints: " , waypoints)
print("waypoints flattened: " , waypoints.flatten())
waypoint_const_violation = path_gen.evaluate_waypoint_constraint(waypoints, control_points_2)
print("waypoint_const_violation: " , waypoint_const_violation)

print(" ")
print("PATH 2 DATA:")
print("computation time path 2: " , end_time_2 - start_time_2)
print("max incline: " , np.max(abs(incline_data_2)))
print("max curvature" , np.max(curvature_data_2))
print(" ")
plt.figure()
ax = plt.axes(projection='3d')
ax.plot(spline_data_1[0,:], spline_data_1[1,:],spline_data_1[2,:], color = "b")
# ax.scatter(control_points[0,:], control_points[1,:],control_points[2,:])
ax.plot(spline_data_2[0,:], spline_data_2[1,:],spline_data_2[2,:], color = "r")
ax.set_xlabel('$X$')
ax.set_ylabel('$Y$')
ax.set_zlabel('$Z$')
plot3D_waypoints(waypoint_data_1, ax)
plot3D_waypoints(waypoint_data_2, ax)
set_axes_equal(ax,dimension)
plt.show()


plt.figure()
plt.plot(spline_data_1[0,:], spline_data_1[1,:], color = "b")
plt.plot(spline_data_2[0,:], spline_data_2[1,:], color = "r")
plt.show()

# plt.figure()
# plt.plot(time_data_1, acceleration_data_1[0,:], color = "b")
# plt.plot(time_data_1, acceleration_data_1[1,:], color = "b")
# plt.plot(time_data_1, acceleration_data_1[2,:], color = "b")
# plt.plot(time_data_2, acceleration_data_2[0,:], color = "r")
# plt.plot(time_data_2, acceleration_data_2[1,:], color = "r")
# plt.plot(time_data_2, acceleration_data_2[2,:], color = "r")
# plt.show()

plt.figure()
plt.plot(time_data_1, curvature_data_1, color = "b")
plt.plot(time_data_2, curvature_data_2, color = "r")
plt.show()