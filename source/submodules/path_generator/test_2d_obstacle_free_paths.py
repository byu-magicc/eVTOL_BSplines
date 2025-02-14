import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from bsplinegenerator.bsplines import BsplineEvaluation
# from path_generation.path_generator import PathGenerator
from path_generation.path_generator import PathGenerator
from path_generation.safe_flight_corridor import SFC, SFC_Data, plot_sfcs, get2DRotationAndTranslationFromPoints
from path_generation.obstacle import Obstacle, plot_2D_obstacles
from path_generation.waypoint_data import plot2D_waypoints
from path_generation.waypoint_data import Waypoint, WaypointData
from path_generation.path_plotter import set_axes_equal
import time


# waypoints case 1
waypoint_a = Waypoint(location=np.array([[0],[0]]), velocity = np.array([[1],[0]]))
waypoint_b = Waypoint(location=np.array([[20],[0]]), velocity = np.array([[1],[0]]))
waypoint_data_1 = WaypointData(start_waypoint=waypoint_a, end_waypoint=waypoint_b)

# waypoints case 2
waypoint_a = Waypoint(location=np.array([[0],[0]]), velocity = np.array([[1],[0]]))
waypoint_b = Waypoint(location=np.array([[20],[-10]]), velocity = np.array([[0],[-1]]))
waypoint_data_2 = WaypointData(start_waypoint=waypoint_a, end_waypoint=waypoint_b)

# waypoints case 3
waypoint_a = Waypoint(location=np.array([[0],[0]]), velocity = np.array([[0],[1]]))
waypoint_b = Waypoint(location=np.array([[20],[0]]), velocity = np.array([[0],[1]]))
waypoint_data_3 = WaypointData(start_waypoint=waypoint_a, end_waypoint=waypoint_b)


waypoint_data_list = [waypoint_data_1, waypoint_data_2, waypoint_data_3]


#obstacle_list 1a
obstacles_1a = None

#obstacle_list 1b
obstacle_1 = Obstacle(center=np.array([[7],[-1]]), radius=1.5)
obstacles_1b = [obstacle_1]

#obstacle_list 1c
obstacle_2 = Obstacle(center=np.array([[15],[1.5]]), radius=2)
obstacles_1c = [obstacle_1,obstacle_2]

#obstacle_list 2a
obstacles_2a = None

#obstacle_list 2b
obstacle_1 = Obstacle(center=np.array([[7],[-2]]), radius=1.5)
obstacles_2b = [obstacle_1]

#obstacle_list 2c
obstacle_2 = Obstacle(center=np.array([[15],[-8]]), radius=2)
obstacles_2c = [obstacle_1,obstacle_2]

#obstacle_list 3a
obstacles_3a = None

#obstacle_list 3b
obstacle_1 = Obstacle(center=np.array([[10],[0]]), radius=1.5)
obstacles_3b = [obstacle_1]

#obstacle_list 3c
obstacle_2 = Obstacle(center=np.array([[5],[2]]), radius=1)
obstacles_3c = [obstacle_1,obstacle_2]


obstacles_list_1 = [obstacles_1a, obstacles_1b, obstacles_1c]
obstacles_list_2 = [obstacles_2a, obstacles_2b, obstacles_2c]
obstacles_list_3 = [obstacles_3a, obstacles_3b, obstacles_3c]
obstacles_list_all = [obstacles_list_1, obstacles_list_2, obstacles_list_3]


order = 3
dimension = 2
num_intervals_free_space = 5
# objective_function_type="minimal_acceleration_path"
objective_function_type="minimal_distance_path"
max_curvature=0.5

fig, ax = plt.subplots(3,3)
for i in range(len(waypoint_data_list)):
    waypoint_data = waypoint_data_list[i]
    obstacles_list = obstacles_list_all[i]
    for j in range(len(obstacles_list_all[i])):
        if j == 2: num_intervals_free_space = 7 
        else: num_intervals_free_space = 5 
        path_gen = PathGenerator(dimension,num_intervals_free_space=num_intervals_free_space)
        obstacles = obstacles_list[j]
        start_time = time.time()
        control_points, status = path_gen.generate_path(waypoint_data=waypoint_data, max_curvature=max_curvature,
            max_incline=None, sfc_data=None, obstacles=obstacles,objective_function_type=objective_function_type)
        end_time = time.time()
        num_cont_pts = np.shape(control_points)[1]
        eval_time = end_time - start_time
        # print("control_points: " , control_points)
        print("computation time:" , end_time - start_time)
        spline_start_time = 0
        scale_factor = 1
        bspline = BsplineEvaluation(control_points, order, spline_start_time, scale_factor, False)
        number_data_points = 10000
        spline_data, time_data = bspline.get_spline_data(number_data_points)
        path_length = round(bspline.get_arc_length(number_data_points),2)
        waypoints = waypoint_data.get_waypoint_locations()
        ax[i,j].plot(spline_data[0,:], spline_data[1,:])
        ax[i,j].set_xlabel("x (m) \n \n evaluation time: " + str(np.round(eval_time,2)) + 
                           "\n path length: " + str(path_length) +
                           "\n num ctrl pts: " + str(num_cont_pts))
        ax[i,j].set_aspect('equal')
        plot2D_waypoints(waypoint_data, ax[i,j],arrow_scale = 2)
        if obstacles is not None:
            plot_2D_obstacles(obstacles, ax[i,j])

# curvature_data, time_data = bspline.get_spline_curvature_data(1000)
# acceleration_data, time_data = bspline.get_spline_derivative_data(1000,2)
# minvo_cps = bspline.get_minvo_control_points()
ax[0,0].set_title("Zero Obstacles")
ax[0,1].set_title("One Obstacle")
ax[0,2].set_title("Two Obstacles")
ax[0,0].set_ylabel("Case 1 \n y (m)")
ax[1,0].set_ylabel("Case 2 \n y (m)")
ax[2,0].set_ylabel("Case 3 \n y (m)")
ax[0,0].set_ylim([-1, 1])
plt.show()

# plt.figure()
# plt.plot(time_data, curvature_data)
# plt.show()

# print("start curvature: " , curvature_data[0])
# print("end curvature: " , curvature_data[-1])


