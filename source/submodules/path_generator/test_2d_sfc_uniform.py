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


#sfc constraints case 
point_1 = np.array([[0],[16]])
point_2 = np.array([[0],[6]])
point_3 = np.array([[7],[6]])
point_4 = np.array([[7],[-1]])
point_5 = np.array([[18],[-1]])
point_6 = np.array([[14],[12]])
point_sequence = np.concatenate((point_1, point_2, point_3, point_4, point_5, point_6),axis=1)
R1, T1, min_len_1 = get2DRotationAndTranslationFromPoints(point_1, point_2)
R2, T2, min_len_2 = get2DRotationAndTranslationFromPoints(point_2, point_3)
R3, T3, min_len_3 = get2DRotationAndTranslationFromPoints(point_3, point_4)
R4, T4, min_len_4 = get2DRotationAndTranslationFromPoints(point_4, point_5)
R5, T5, min_len_5 = get2DRotationAndTranslationFromPoints(point_5, point_6)
sfc_1 = SFC(np.array([[min_len_1+4],[3]]), T1, R1)
sfc_2 = SFC(np.array([[min_len_2+3],[4]]), T2, R2)
sfc_3 = SFC(np.array([[min_len_3+4],[3]]), T3, R3)
sfc_4 = SFC(np.array([[min_len_4+3],[3]]), T4, R4)
sfc_5 = SFC(np.array([[min_len_5+3],[3]]), T5, R5)
sfcs = [sfc_1, sfc_2, sfc_3, sfc_4, sfc_5]
sfc_data = SFC_Data(sfcs,point_sequence)
obstacle = Obstacle(center=np.array([[14],[0]]), radius=1.5)
obstacles = [obstacle]
obstacles = None
waypoint_a = Waypoint(location=point_1, velocity = (point_2-point_1)/np.linalg.norm((point_2-point_1)))
waypoint_b = Waypoint(location=point_6, velocity = (point_6-point_5)/np.linalg.norm((point_6-point_5)))
waypoint_data = WaypointData(start_waypoint=waypoint_a, end_waypoint=waypoint_b)


order = 3
dimension = 2
num_intervals_free_space = 5
# objective_function_type="minimal_acceleration_path"
objective_function_type = "minimal_distance_path"
max_curvature=0.5

fig, ax = plt.subplots()
path_gen = PathGenerator(dimension,num_intervals_free_space=num_intervals_free_space)
start_time = time.time()
control_points, status = path_gen.generate_path(waypoint_data=waypoint_data, max_curvature=max_curvature,
        max_incline=None, sfc_data=sfc_data, obstacles=obstacles,objective_function_type=objective_function_type)
end_time = time.time()
eval_time = end_time - start_time
num_cont_pts = np.shape(control_points)[1]
print("computation time:" , end_time - start_time)
spline_start_time = 0
scale_factor = 1
bspline = BsplineEvaluation(control_points, order, spline_start_time, scale_factor, False)
number_data_points = 10000
spline_data, time_data = bspline.get_spline_data(number_data_points)
path_length = bspline.get_arc_length(number_data_points)
waypoints = waypoint_data.get_waypoint_locations()
ax.plot(spline_data[0,:], spline_data[1,:])
ax.set_xlabel("x (m) \n evaluation time: " + str(np.round(eval_time,2)) + "\n num ctrl pts: " + str(num_cont_pts))
ax.set_aspect('equal')
plot2D_waypoints(waypoint_data, ax, arrow_scale = 2)

alpha = 1
plot_sfcs(sfcs, ax, alpha=alpha)
if obstacles is not None:
    plot_2D_obstacles(obstacles, ax)


plt.show()




