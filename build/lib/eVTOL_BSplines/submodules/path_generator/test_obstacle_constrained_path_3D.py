import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from bsplinegenerator.bsplines import BsplineEvaluation
from path_generation.path_generator import PathGenerator
from path_generation.obstacle import Obstacle, plot_3D_obstacles
from path_generation.path_plotter import set_axes_equal
from path_generation.waypoint_data import Waypoint, WaypointData
import time

point_1 = np.array([[3],[4],[0]])
point_2 = np.array([[7],[10],[13]])
points = np.concatenate((point_1,point_2),1)
dimension = np.shape(point_1)[0]

waypoint_1 = Waypoint(location=point_1)
waypoint_2 = Waypoint(location=point_2)
# waypoint_1.velocity = point_2- point_1
# waypoint_2.velocity = np.array([[0],[1],[0]])
waypoint_data = WaypointData(start_waypoint=waypoint_1,end_waypoint=waypoint_2)

max_curvature = 0.5
max_incline = 2
obstacles = [Obstacle(np.array([[4.5],[5.5],[3]]), 1.0),
             Obstacle(np.array([[4],[8],[9.5]]), 1.5)]

order = 3

path_gen = PathGenerator(dimension)
start_time = time.time()
control_points, status = path_gen.generate_path(waypoint_data=waypoint_data, max_curvature=max_curvature,
    max_incline=max_incline, sfc_data=None, obstacles=obstacles)
end_time = time.time()
print("control_points: " , control_points)
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
acceleration_data, time_data = bspline.get_spline_derivative_data(number_data_points,2)
# acceleration_data, time_data = bspline.get_spline_derivative_data(1000,2)
minvo_cps = bspline.get_minvo_control_points()
waypoints = waypoint_data.get_waypoint_locations()

print("max incline" , np.max(incline_data))
print("max curvature" , np.max(curvature_data))
print("computation time: " , end_time - start_time)

cp_accel = (control_points[:,0] - 2*control_points[:,1] + control_points[:,2])
# print("cp_accel:  " , cp_accel)
# print("start accel: " , acceleration_data[:,0])
# print("end accel: " , acceleration_data[:,-1])

# print("sfcs: " , sfcs)
plt.figure()
ax = plt.axes(projection='3d')
if obstacles != None:
    plot_3D_obstacles(obstacles, ax)
ax.plot(spline_data[0,:], spline_data[1,:],spline_data[2,:])
ax.scatter(waypoints[0,:],waypoints[1,:],waypoints[2,:])
ax.set_xlabel('$X$')
ax.set_ylabel('$Y$')
ax.set_zlabel('$Z$')
set_axes_equal(ax, dimension)
ax.set_xlabel('$X$')
ax.set_ylabel('$Y$')
ax.set_zlabel('$Z$')
# plt.scatter(minvo_cps[0,:],minvo_cps[1,:])
plt.show()



# plt.figure()
# plt.plot(time_data, incline_data)
# plt.show()

# plt.figure()
# plt.plot(time_data, curvature_data)
# plt.show()


