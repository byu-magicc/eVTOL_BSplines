import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from bsplinegenerator.bsplines import BsplineEvaluation
from trajectory_generation.trajectory_generator import TrajectoryGenerator
from trajectory_generation.constraint_data_structures.safe_flight_corridor import SFC_Data, SFC, plot_sfcs, get2DRotationAndTranslationFromPoints
from trajectory_generation.path_plotter import set_axes_equal
from trajectory_generation.constraint_data_structures.waypoint_data import Waypoint, WaypointData, plot2D_waypoints
from trajectory_generation.constraint_data_structures.dynamic_bounds import DerivativeBounds, TurningBound
from trajectory_generation.constraint_data_structures.obstacle import Obstacle, plot_2D_obstacles
from trajectory_generation.constraint_data_structures.constraints_container import ConstraintsContainer
import time

#note incline constraints work much better when have a start and an end direction.

dimension = 2
# max_curvature = 1
order = 3
# traj_objective_type = "minimal_acceleration_path"
# traj_objective_type = "minimal_velocity_path"
# traj_objective_type = "minimal_distance_path"
traj_objective_type = "minimal_time_path"

waypoint_1 = Waypoint(location=np.array([[-5],[0]]),velocity=np.array([[0],[15]]))
waypoint_2 = Waypoint(location=np.array([[5],[0]]),velocity=np.array([[0],[10]]))

point_1 = np.array([[-5],[0]])
point_2 = np.array([[0],[5]])
point_3 = np.array([[0],[-5]])
point_4 = np.array([[5],[0]])
point_sequence = np.concatenate((point_1,point_2,point_3,point_4),1)
dimension = np.shape(point_1)[0]
R1, T1, min_len_1 = get2DRotationAndTranslationFromPoints(point_1, point_2)
R2, T2, min_len_2 = get2DRotationAndTranslationFromPoints(point_2, point_3)
R3, T3, min_len_3 = get2DRotationAndTranslationFromPoints(point_3, point_4)
sfc_1 = SFC(np.array([[min_len_1+3],[2]]), T1, R1)
sfc_2 = SFC(np.array([[min_len_2 + 2],[3]]), T2, R2)
sfc_3 = SFC(np.array([[min_len_3+3],[2]]), T3, R3)
sfcs = (sfc_1, sfc_2, sfc_3)
min_intervals_per_corridor = 1
sfc_data = SFC_Data(sfcs, point_sequence,min_intervals_per_corridor,intervals_per_corridor=np.array([1,1,1]))

# sfc_data = None
# obstacles = [Obstacle(center=np.array([[5.5],[7]]), radius=1)]
obstacles = None
obstacle_list = None

# max_turning_bound = 1 #angular rate
# turning_bound = TurningBound(max_turning_bound,"angular_rate")

# max_turning_bound = 0.5 #cent accel
# turning_bound = TurningBound(max_turning_bound,"centripetal_acceleration")

# max_turning_bound = 1 #curv
# turning_bound = TurningBound(max_turning_bound,"curvature")

turning_bound = None

max_velocity = 30
# min_velocity = 4
# max_velocity = None
max_acceleration = 100
# max_acceleration = 30
# max_acceleration = None
derivative_bounds = DerivativeBounds(max_velocity, max_acceleration)
# derivative_bounds = DerivativeBounds(max_velocity = max_velocity, 
                                    #  min_velocity=min_velocity)
# derivative_bounds = None


waypoint_sequence = (waypoint_1, waypoint_2)
waypoint_data = WaypointData(waypoint_sequence)
traj_gen = TrajectoryGenerator(dimension)
start_time_1 = time.time()


constraints_container = ConstraintsContainer(waypoint_constraints = waypoint_data, derivative_constraints=derivative_bounds,
    turning_constraint=turning_bound, sfc_constraints=sfc_data, obstacle_constraints=obstacle_list)

control_points, scale_factor, is_violation = traj_gen.generate_trajectory(constraints_container, traj_objective_type, num_intervals_free_space=10)
end_time_1 = time.time()
spline_start_time_1 = 0
bspline = BsplineEvaluation(control_points, order, spline_start_time_1, scale_factor, False)
end_time_spline = bspline.get_end_time()

# bspline.plot_spline(1000)

## spline 1 data
number_data_points = 10000
spline_data, time_data = bspline.get_spline_data(number_data_points)
curvature_data, time_data = bspline.get_spline_curvature_data(number_data_points)
velocity_data, time_data = bspline.get_derivative_magnitude_data(number_data_points,1)
acceleration_data, time_data = bspline.get_derivative_magnitude_data(number_data_points,2)
angular_rate_data, time_data = bspline.get_angular_rate_data(number_data_points)
centripetal_acceleration_data, time_data = bspline.get_centripetal_acceleration_data(number_data_points)
path_length = bspline.get_arc_length(number_data_points)
start_velocity = bspline.get_derivative_at_time_t(0,1)
start_acceleration = bspline.get_derivative_at_time_t(0,2)
print("path_length: " , path_length)
print("computation time: " , end_time_1 - start_time_1)

velocity_matrix, time_data = bspline.get_spline_derivative_data(number_data_points,1)
acceleration_matrix, time_data = bspline.get_spline_derivative_data(number_data_points,2)

plt.figure()
ax = plt.axes()
# ax.scatter(control_points[0,:], control_points[1,:], color="tab:orange")
ax.plot(spline_data[0,:], spline_data[1,:], color = "tab:blue", label="optimized path")
plot2D_waypoints(waypoint_data, ax)
plot_sfcs(sfcs, ax)
set_axes_equal(ax,dimension)
plt.xlim(np.min(spline_data[0,:])-2, np.max(spline_data[0,:])+2)
plt.ylim(np.min(spline_data[1,:])-2, np.max(spline_data[1,:])+2)
plt.title("B-spline Trajectory")
plt.xlabel("x position (m)")
plt.ylabel("y position (m)")
plt.legend()
plt.show()

# if max_velocity is not None:
#     plt.figure()
#     plt.plot(time_data, velocity_data,color = "b")
#     plt.plot(time_data, max_velocity + velocity_data*0)
#     plt.plot(time_data, min_velocity + velocity_data*0)
#     plt.title("velocity")
#     plt.show()

# if max_acceleration is not None:
#     plt.figure()
#     plt.plot(time_data, acceleration_data,color = "b")
#     plt.plot(time_data, max_acceleration + acceleration_data*0)
#     plt.title("acceleration")
#     plt.show()

# if turning_bound is not None:
#     turn_data = []
#     if turning_bound.bound_type == "angular_rate":
#         turn_data = angular_rate_data
#     elif turning_bound.bound_type == "curvature":
#         turn_data = curvature_data
#     elif turning_bound.bound_type == "centripetal_acceleration":
#         turn_data = centripetal_acceleration_data
#     turn_title = turning_bound.bound_type
#     plt.figure()
#     plt.plot(time_data, turn_data,color = "b")
#     # plt.plot(time_data, acceleration_data,color = "g")
#     plt.plot(time_data, max_turning_bound + turn_data*0)
#     plt.title(turn_title)
#     plt.show()


# create 2 subplots
fig, ax = plt.subplots(nrows=2, ncols=1)
fig.suptitle("Derivative Data")
#velocity data
ax[0].set_ylabel("velocity (m/s)")
ax[0].plot(time_data, velocity_data, label="velocity magnitude")
ax[0].scatter([time_data[0], time_data[-1]], [velocity_data[0],velocity_data[-1]],
              facecolors='none', edgecolors='r', label ="terminal velocity constraints")
ax[0].plot(time_data, velocity_data*0+max_velocity, color='r', label="max velocity")
ax[0].legend()
# facecolors='none', edgecolors='r'
#acceleration data
ax[1].set_ylabel("acceleration (m/s^2)")
ax[1].plot(time_data, acceleration_data, label="acceleration magnitude")
ax[1].scatter([time_data[0], time_data[-1]], [acceleration_data[0],acceleration_data[-1]],
              facecolors='none', edgecolors='r', label= "terminal acceleration constraints")
ax[1].plot(time_data, acceleration_data*0+max_acceleration, color='r', label="max acceleration")
ax[1].legend()

plt.show()
