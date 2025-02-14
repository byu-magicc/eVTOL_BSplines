import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from bsplinegenerator.bsplines import BsplineEvaluation
# from path_generation.path_generator import PathGenerator
from path_generation.path_generator import PathGenerator
from path_generation.safe_flight_corridor import SFC_Data, get3DRotationAndTranslationFromPoints
from path_generation.path_plotter import set_axes_equal
from path_generation.waypoint_data import Waypoint, WaypointData, plot3D_waypoints
from path_generation.obstacle import Obstacle, plot_3D_obstacles, plot_3D_obstacle
import time

def set_velocity_scale(waypoint_a: Waypoint, waypoint_b: Waypoint, custom_scale:float=None):
    distance = np.linalg.norm(waypoint_a.location - waypoint_b.location)
    if custom_scale is None:
        scale = distance/1.1
    else:
        scale = custom_scale
    vel_norm_a = np.linalg.norm(waypoint_a.velocity)
    vel_norm_b = np.linalg.norm(waypoint_b.velocity)
    waypoint_a.velocity = waypoint_a.velocity/vel_norm_a * scale
    waypoint_a.acceleration = waypoint_a.acceleration/(vel_norm_a**2) *scale**2
    waypoint_b.velocity = waypoint_b.velocity/vel_norm_b *scale

##### Initial Settings
dimension = 3
max_curvature = 0.0114362
r_min = 1/max_curvature
max_incline = 0.2
order = 3
scale_factor = 1
path_objective_type = "minimal_velocity_path"
obstacle_type="sphere"
obstacle_radius = 50
obstacle_center = np.array([[250],[250],[150]])
obstacle = Obstacle(center=obstacle_center, radius=obstacle_radius)
obstacles = [obstacle] 
obstacles = None
number_data_points = 10000
num_intervals_free_space = 9
path_gen = PathGenerator(dimension, num_intervals_free_space=num_intervals_free_space)
scale = None
scale = 600

### 1st path
acceleration = np.array([[0],[0],[0]])
waypoint_0 = Waypoint(location = np.array([600,-100,300])[:,None],
                      velocity = np.array([-600,  0, 0])[:,None],
                      acceleration = np.array([0,0,0])[:,None])
waypoint_1 = Waypoint(location = np.array([155.82943463, 483.11508766,  91.58818684])[:,None],
                      velocity = np.array([218.30452428,  98.65421414,   0.])[:,None],
                      acceleration = np.array([0.88817539, -4.86019262, 52.21121555])[:,None])
# set_velocity_scale(waypoint_0, waypoint_1, scale)
waypoint_data = WaypointData(start_waypoint=waypoint_0,end_waypoint=waypoint_1)
start_time = time.time()
control_points, success = path_gen.generate_path(waypoint_data=waypoint_data, max_curvature=max_curvature,
    max_incline=max_incline, sfc_data=None, obstacles=obstacles, objective_function_type=path_objective_type, 
    obstacle_type=obstacle_type)
end_time = time.time()
path_time = end_time - start_time
bspline = BsplineEvaluation(control_points, order, 0, scale_factor, False)
time_spline_1 = bspline.get_end_time()
print("acceleration: ", bspline.get_derivative_at_time_t(time_spline_1,2))



spline_data, time_data = bspline.get_spline_data(number_data_points)
curvature_data, time_data = bspline.get_spline_curvature_data(number_data_points)
incline_data, time_data = bspline.get_spline_slope_data(number_data_points)


print("computation time: " , path_time)
print("max_curvature: ", np.max(curvature_data))
print("max_slope: ", np.max(np.abs(incline_data)))
print("scale: " , np.linalg.norm(waypoint_0.velocity))

np.save("control_points_0", control_points)


plt.figure()
ax = plt.axes(projection='3d')
ax.plot(spline_data[0,:], spline_data[1,:],spline_data[2,:], color = "b")
ax.set_xlabel('$X$')
ax.set_ylabel('$Y$')
ax.set_zlabel('$Z$')
plot_3D_obstacle(obstacle, ax)
plot3D_waypoints(waypoint_data, ax)
set_axes_equal(ax,dimension)
plt.show()

# plt.figure()
# plt.title("Curvature")
# plt.plot(time_data_1, curvature_data_1,color = "b")
# plt.plot(time_data_2, curvature_data_2, color = "r")
# plt.plot(time_data_3, curvature_data_3, color = "g")
# plt.show()


