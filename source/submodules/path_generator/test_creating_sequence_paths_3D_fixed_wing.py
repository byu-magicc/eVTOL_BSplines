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
# obstacles = None
number_data_points = 10000
num_intervals_free_space = 8
path_gen = PathGenerator(dimension, num_intervals_free_space=num_intervals_free_space)

### zeroth path
waypoint_start = Waypoint(location = np.array([600,-100,300])[:,None],
                         velocity = np.array([-600,  0, 0])[:,None],
                         acceleration = np.array([0,0,0])[:,None])
waypoint_0 = Waypoint(location = np.array([155.82943463, 483.11508766,  91.58818684])[:,None],
                      velocity = np.array([218.30452428,  98.65421414,   0.])[:,None],
                      acceleration = np.array([0.88817539, -4.86019262, 52.21121555])[:,None])
waypoint_data_0 = WaypointData(start_waypoint=waypoint_start,end_waypoint=waypoint_0)
start_time_0 = time.time()
control_points_0, success_0 = path_gen.generate_path(waypoint_data=waypoint_data_0, max_curvature=max_curvature,
    max_incline=max_incline, sfc_data=None, obstacles=obstacles, objective_function_type=path_objective_type, 
    obstacle_type=obstacle_type)
end_time_0 = time.time()
np.save("control_points_0", control_points_0)
path_time_0 = end_time_0 - start_time_0
bspline_0 = BsplineEvaluation(control_points_0, order, 0, scale_factor, False)
time_spline_0 = bspline_0.get_end_time()
### 1st path
Waypoint(location = np.array([155.82943463, 483.11508766,  91.58818684])[:,None],
                      velocity = np.array([218.30452428,  98.65421414,   0.])[:,None],
                      acceleration = np.array([0.88817539, -4.86019262, 52.21121555])[:,None])
waypoint_1 = Waypoint(location=np.array([25.6464104,  151.96743085,  41.38830642])[:,None],
                      velocity=np.array([208.20197419, 133.55090968,   0.])[:,None],
                      acceleration=np.array([-365.70954522,  581.71000345,   78.08153607])[:,None])
waypoint_data_1 = WaypointData(start_waypoint=waypoint_0,end_waypoint=waypoint_1)
start_time_1 = time.time()
control_points_1, success_1 = path_gen.generate_path(waypoint_data=waypoint_data_1, max_curvature=max_curvature,
    max_incline=max_incline, sfc_data=None, obstacles=obstacles, objective_function_type=path_objective_type, 
    obstacle_type=obstacle_type)
end_time_1 = time.time()
np.save("control_points_1", control_points_1)
path_time_1 = end_time_1 - start_time_1
bspline_1 = BsplineEvaluation(control_points_1, order, 0, scale_factor, False)
time_spline_1 = bspline_1.get_end_time()


#### 2nd path 
waypoint_1 = Waypoint(location=np.array([25.6464104,  151.96743085,  41.38830642])[:,None],
                      velocity=np.array([208.20197419, 133.55090968,   0.])[:,None],
                      acceleration=np.array([-365.70954522,  581.71000345,   78.08153607])[:,None])
waypoint_2 = Waypoint(location=np.array([334.66865557, 300.21520007, 183.48560709])[:,None],
                      velocity=np.array([129.60095992,  50.08801435,   0.])[:,None],
                      acceleration=np.array([-62.08330637, 197.34112296, -65.54292566])[:,None])
waypoint_data_2 = WaypointData(start_waypoint=waypoint_1, end_waypoint=waypoint_2)
start_time_2 = time.time()
control_points_2, success_2 = path_gen.generate_path(waypoint_data=waypoint_data_2, max_curvature=max_curvature,
    max_incline=max_incline, sfc_data=None, obstacles=obstacles, objective_function_type=path_objective_type, 
    obstacle_type=obstacle_type)
end_time_2 = time.time()
np.save("control_points_2", control_points_2)
path_time_2 = end_time_2 - start_time_2
bspline_2 = BsplineEvaluation(control_points_2, order, 0, scale_factor, False)
time_spline_2 = bspline_2.get_end_time()

##### 3rd path
waypoint_2 = Waypoint(location=np.array([334.66865557, 300.21520007, 183.48560709])[:,None],
                      velocity=np.array([129.60095992,  50.08801435,   0.])[:,None],
                      acceleration=np.array([-62.08330637, 197.34112296, -65.54292566])[:,None])
waypoint_3 = Waypoint(location=np.array([174.25495117, 349.51711852,  59.89958539])[:,None],
                      velocity=np.array([93.07189172, 179.39280427,   0.])[:,None],
                      acceleration=np.array([193.28489124,  17.7495171,   -7.70208549])[:,None])
waypoint_data_3 = WaypointData(start_waypoint=waypoint_2, end_waypoint=waypoint_3)
start_time_3 = time.time()
control_points_3, success_3 = path_gen.generate_path(waypoint_data=waypoint_data_3, max_curvature=max_curvature,
    max_incline=max_incline, sfc_data=None, obstacles=obstacles, objective_function_type=path_objective_type, 
    obstacle_type=obstacle_type)
end_time_3 = time.time()
np.save("control_points_3", control_points_3)
path_time_3 = end_time_3 - start_time_3
bspline_3 = BsplineEvaluation(control_points_3, order, 0, scale_factor, False)
time_spline_3 = bspline_3.get_end_time()


##### 4th path
waypoint_3 = Waypoint(location=np.array([174.25495117, 349.51711852,  59.89958539])[:,None],
                      velocity=np.array([93.07189172, 179.39280427,   0.])[:,None],
                      acceleration=np.array([193.28489124,  17.7495171,   -7.70208549])[:,None])
waypoint_4 = Waypoint(location=np.array([441.52554136, 491.61093362,  43.31258304])[:,None],
                      velocity=np.array([91.44251366, 83.33792852,  0.])[:,None],
                      acceleration=np.array([-38.14621985,  42.75159813,   9.02912899])[:,None])
waypoint_data_4 = WaypointData(start_waypoint=waypoint_3, end_waypoint=waypoint_4)
start_time_4 = time.time()
control_points_4, success_4 = path_gen.generate_path(waypoint_data=waypoint_data_4, max_curvature=max_curvature,
    max_incline=max_incline, sfc_data=None, obstacles=obstacles, objective_function_type=path_objective_type, 
    obstacle_type=obstacle_type)
end_time_4 = time.time()
np.save("control_points_4", control_points_4)
path_time_4 = end_time_4 - start_time_4
bspline_4 = BsplineEvaluation(control_points_4, order, 0, scale_factor, False)
time_spline_4 = bspline_4.get_end_time()


##### 5th path
waypoint_4 = Waypoint(location=np.array([441.52554136, 491.61093362,  43.31258304])[:,None],
                      velocity=np.array([91.44251366, 83.33792852,  0.])[:,None],
                      acceleration=np.array([-38.14621985,  42.75159813,   9.02912899])[:,None])
waypoint_5 = Waypoint(location=np.array([467.4974241,  307.97318145,  49.88985147])[:,None],
                      velocity=np.array([78.30783168, 95.78531794,  0.])[:,None])
waypoint_data_5 = WaypointData(start_waypoint=waypoint_4, end_waypoint=waypoint_5)
start_time_5 = time.time()
control_points_5, success_5 = path_gen.generate_path(waypoint_data=waypoint_data_5, max_curvature=max_curvature,
    max_incline=max_incline, sfc_data=None, obstacles=obstacles, objective_function_type=path_objective_type, 
    obstacle_type=obstacle_type)
end_time_5 = time.time()
np.save("control_points_5", control_points_5)
path_time_5 = end_time_5 - start_time_5
bspline_5 = BsplineEvaluation(control_points_5, order, 0, scale_factor, False)
time_spline_5 = bspline_5.get_end_time()
spline_data_5, time_data_5 = bspline_5.get_spline_data(number_data_points)

waypoints = np.concatenate((waypoint_0.location, waypoint_1.location,waypoint_2.location, waypoint_3.location, waypoint_4.location, waypoint_5.location),1)
np.save("waypoints" , waypoints)

spline_data_0, time_data_0 = bspline_0.get_spline_data(number_data_points)
spline_data_1, time_data_1 = bspline_1.get_spline_data(number_data_points)
spline_data_2, time_data_2 = bspline_2.get_spline_data(number_data_points)
spline_data_3, time_data_3 = bspline_3.get_spline_data(number_data_points)
spline_data_4, time_data_4 = bspline_4.get_spline_data(number_data_points)
spline_data_5, time_data_5 = bspline_5.get_spline_data(number_data_points)

curvature_data_0, time_data_0 = bspline_0.get_spline_curvature_data(number_data_points)
curvature_data_1, time_data_1 = bspline_1.get_spline_curvature_data(number_data_points)
curvature_data_2, time_data_2 = bspline_2.get_spline_curvature_data(number_data_points)
curvature_data_3, time_data_3 = bspline_3.get_spline_curvature_data(number_data_points)
curvature_data_4, time_data_4 = bspline_4.get_spline_curvature_data(number_data_points)
curvature_data_5, time_data_5 = bspline_5.get_spline_curvature_data(number_data_points)

incline_data_0, time_data_0 = bspline_0.get_spline_slope_data(number_data_points)
incline_data_1, time_data_1 = bspline_1.get_spline_slope_data(number_data_points)
incline_data_2, time_data_2 = bspline_2.get_spline_slope_data(number_data_points)
incline_data_3, time_data_3 = bspline_3.get_spline_slope_data(number_data_points)
incline_data_4, time_data_4 = bspline_4.get_spline_slope_data(number_data_points)
incline_data_5, time_data_5 = bspline_5.get_spline_slope_data(number_data_points)

waypoint_start.print_waypoint(-1)
waypoint_0.print_waypoint(0)
waypoint_1.print_waypoint(1)
waypoint_2.print_waypoint(2)
waypoint_3.print_waypoint(3)
waypoint_4.print_waypoint(4)
waypoint_5.print_waypoint(5)

print("computation time 0: " , path_time_0)
if (max_curvature < np.max(curvature_data_0)):
    print(" max curvature 0: ", np.max(curvature_data_0))
if (max_incline < np.max(np.abs(incline_data_0))):
    print(" max_slope 0: " , np.max(np.abs(incline_data_0)))

print("computation time 1: " , path_time_1)
if (max_curvature < np.max(curvature_data_1)):
    print(" max curvature 1: ", np.max(curvature_data_1))
if (max_incline < np.max(np.abs(incline_data_1))):
    print(" max_slope 1: " , np.max(np.abs(incline_data_1)))

print("computation time 2: " , path_time_2)
if (max_curvature < np.max(curvature_data_2)):
    print(" max curvature 2: ", np.max(curvature_data_2))
if (max_incline < np.max(np.abs(incline_data_2))):
    print(" max_slope 2: " , np.max(np.abs(incline_data_2)))

print("computation time 3: " , path_time_3)
if (max_curvature < np.max(curvature_data_3)):
    print(" max curvature 3: ", np.max(curvature_data_3))
if (max_incline < np.max(np.abs(incline_data_3))):
    print(" max_slope 3: " , np.max(np.abs(incline_data_3)))

print("computation time 4: " , path_time_4)
if (max_curvature < np.max(curvature_data_4)):
    print(" max curvature 4: ", np.max(curvature_data_4))
if (max_incline < np.max(np.abs(incline_data_4))):
    print(" max_slope 4: " , np.max(np.abs(incline_data_4)))

print("computation time 5: " , path_time_5)
if (max_curvature < np.max(curvature_data_4)):
    print(" max curvature 5: ", np.max(curvature_data_5))
if (max_incline < np.max(np.abs(incline_data_5))):
    print(" max_slope 5: " , np.max(np.abs(incline_data_5)))

label_0 = "path 0: " + str(np.round(path_time_0,2))  + " sec"
label_1 = "path 1: " + str(np.round(path_time_1,2))  + " sec"
label_2 = "path 2: " + str(np.round(path_time_2,2))  + " sec"
label_3 = "path 3: " + str(np.round(path_time_3,2))  + " sec"
label_4 = "path 4: " + str(np.round(path_time_4,2))  + " sec"
label_5 = "path 5: " + str(np.round(path_time_5,2))  + " sec"
                           
                           
plt.figure()
ax = plt.axes(projection='3d')
ax.plot(spline_data_0[0,:], spline_data_0[1,:],spline_data_0[2,:], color = "m", label=label_0)
ax.plot(spline_data_1[0,:], spline_data_1[1,:],spline_data_1[2,:], color = "b", label=label_1)
ax.plot(spline_data_2[0,:], spline_data_2[1,:],spline_data_2[2,:], color = "r", label=label_2)
ax.plot(spline_data_3[0,:], spline_data_3[1,:],spline_data_3[2,:], color = "g", label=label_3)
ax.plot(spline_data_4[0,:], spline_data_4[1,:],spline_data_4[2,:], color = "tab:olive", label=label_4)
ax.plot(spline_data_5[0,:], spline_data_5[1,:],spline_data_5[2,:], color = "tab:orange", label=label_5)
ax.set_xlabel('X(m)')
ax.set_ylabel(' Y(m)')
ax.set_zlabel(' Z(m)')
arrow_length = 100
plot_3D_obstacle(obstacle, ax)
plot3D_waypoints(waypoint_data_0, ax, arrow_length = arrow_length)
plot3D_waypoints(waypoint_data_1, ax, arrow_length = arrow_length)
# plot3D_waypoints(waypoint_data_2, ax)
plot3D_waypoints(waypoint_data_3, ax, arrow_length = arrow_length)
# plot3D_waypoints(waypoint_data_4, ax)
plot3D_waypoints(waypoint_data_5, ax, arrow_length = arrow_length)
set_axes_equal(ax,dimension)
# ax.legend()
plt.show()

plt.figure()
plt.title("Curvature")
plt.plot(time_data_1, curvature_data_1,color = "b")
plt.plot(time_data_2+time_data_1[-1], curvature_data_2, color = "r")
plt.plot(time_data_3+time_data_2[-1]+time_data_1[-1], curvature_data_3, color = "g")
plt.plot(time_data_4+time_data_2[-1]+time_data_1[-1]+time_data_3[-1], curvature_data_4, color = "tab:olive")
plt.plot(time_data_5+time_data_2[-1]+time_data_1[-1]+time_data_3[-1]+time_data_4[-1], curvature_data_5, color = "tab:orange")
plt.plot(np.array([0,time_data_5[-1]+time_data_2[-1]+time_data_1[-1]+time_data_3[-1]+time_data_4[-1]]), np.array([1,1])*max_curvature, color = 'k')
plt.show()

plt.figure()
plt.title("Slope")
plt.plot(time_data_1, np.abs(incline_data_1),color = "b")
plt.plot(time_data_2+time_data_1[-1], np.abs(incline_data_2), color = "r")
plt.plot(time_data_3+time_data_2[-1]+time_data_1[-1], np.abs(incline_data_3), color = "g")
plt.plot(time_data_4+time_data_2[-1]+time_data_1[-1]+time_data_3[-1], np.abs(incline_data_4), color = "tab:olive")
plt.plot(time_data_5+time_data_2[-1]+time_data_1[-1]+time_data_3[-1]+time_data_4[-1], np.abs(incline_data_5), color = "tab:orange")
plt.plot(np.array([0,time_data_5[-1]+time_data_2[-1]+time_data_1[-1]+time_data_3[-1]+time_data_4[-1]]), np.array([1,1])*max_incline, color = 'k')
plt.show()
