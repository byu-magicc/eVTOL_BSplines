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


#obstacle_list 1a
sfc_data_1a = None

#obstacle_list 1b
point_1 = np.array([[0],[0]])
point_2 = np.array([[6],[0]])
point_3 = np.array([[6],[7]])
point_4 = np.array([[-1],[7]])
point_sequence = np.concatenate((point_1,point_2,point_3,point_4),axis=1)
R1, T1, min_len_1 = get2DRotationAndTranslationFromPoints(point_1, point_2)
R2, T2, min_len_2 = get2DRotationAndTranslationFromPoints(point_2, point_3)
R3, T3, min_len_3 = get2DRotationAndTranslationFromPoints(point_3, point_4)
sfc_1 = SFC(np.array([[min_len_1+4],[3]]), T1, R1)
sfc_2 = SFC(np.array([[min_len_2+3],[4]]), T2, R2)
sfc_3 = SFC(np.array([[min_len_3+4],[3]]), T3, R3)
sfcs = [sfc_1, sfc_2, sfc_3]
sfc_data_1b = SFC_Data(sfcs,point_sequence)

# waypoints case 1
waypoint_a = Waypoint(location=point_1, velocity = (point_2-point_1)/np.linalg.norm((point_2-point_1)))
waypoint_b = Waypoint(location=point_4, velocity = (point_4-point_3)/np.linalg.norm((point_4-point_3)))
waypoint_data_1 = WaypointData(start_waypoint=waypoint_a, end_waypoint=waypoint_b)


#obstacle_list 2a
sfc_data_2a = None

#obstacle_list 2b
point_1 = np.array([[3],[4]])
point_2 = np.array([[7],[10]])
point_3 = np.array([[14],[7]])
point_4 = np.array([[20],[14]])
point_5 = np.array([[25],[5]])
point_sequence = np.concatenate((point_1,point_2,point_3,point_4,point_5),axis=1)
R1, T1, min_len_1 = get2DRotationAndTranslationFromPoints(point_1, point_2)
R2, T2, min_len_2 = get2DRotationAndTranslationFromPoints(point_2, point_3)
R3, T3, min_len_3 = get2DRotationAndTranslationFromPoints(point_3, point_4)
R4, T4, min_len_4 = get2DRotationAndTranslationFromPoints(point_4, point_5)
sfc_1 = SFC(np.array([[min_len_1+5],[3]]), T1, R1)
sfc_2 = SFC(np.array([[min_len_2+3],[4]]), T2, R2)
sfc_3 = SFC(np.array([[min_len_3+5],[3]]), T3, R3)
sfc_4 = SFC(np.array([[min_len_4+5],[4]]), T4, R4)
sfcs = [sfc_1, sfc_2, sfc_3,sfc_4]
sfc_data_2b = SFC_Data(sfcs,point_sequence)

# waypoints case 2
waypoint_a = Waypoint(location=point_1, velocity = (point_2-point_1)/np.linalg.norm((point_2-point_1)))
waypoint_b = Waypoint(location=point_5, velocity = (point_5-point_4)/np.linalg.norm((point_5-point_4)))
waypoint_data_2 = WaypointData(start_waypoint=waypoint_a, end_waypoint=waypoint_b)

#obstacle_list 2a
sfc_data_3a = None

#obstacle_list 1b
point_1 = np.array([[0],[0]])
point_2 = np.array([[0],[6]])
point_3 = np.array([[7],[6]])
point_4 = np.array([[7],[-1]])
point_5 = np.array([[18],[-1]])
point_6 = np.array([[14],[12]])
point_sequence = np.concatenate((point_1,point_2,point_3,point_4,point_5,point_6),axis=1)
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
sfc_data_3b = SFC_Data(sfcs,point_sequence)
# waypoints case 3
waypoint_a = Waypoint(location=point_1, velocity = (point_2-point_1)/np.linalg.norm((point_2-point_1)))
waypoint_b = Waypoint(location=point_6, velocity = (point_6-point_5)/np.linalg.norm((point_6-point_5)))
waypoint_data_3 = WaypointData(start_waypoint=waypoint_a, end_waypoint=waypoint_b)

waypoint_data_list = [waypoint_data_1, waypoint_data_2, waypoint_data_3]

sfc_list_1 = [sfc_data_1a, sfc_data_1b]
sfc_list_2 = [sfc_data_2a, sfc_data_2b]
sfc_list_3 = [sfc_data_3a, sfc_data_3b]
sfc_list_all = [sfc_list_1, sfc_list_2, sfc_list_3]

order = 3
dimension = 2
num_intervals_free_space = 5
# objective_function_type="minimal_acceleration_path"
objective_function_type = "minimal_distance_path"
max_curvature=0.5

fig, ax = plt.subplots(3,2)
for i in range(len(waypoint_data_list)):
    waypoint_data = waypoint_data_list[i]
    sfc_list = sfc_list_all[i]
    for j in range(len(sfc_list_all[i])):
        path_gen = PathGenerator(dimension,num_intervals_free_space=num_intervals_free_space)
        sfc_data = sfc_list[j]
        start_time = time.time()
        control_points, status = path_gen.generate_path(waypoint_data=waypoint_data, max_curvature=max_curvature,
            max_incline=None, sfc_data=sfc_data, obstacles=None,objective_function_type=objective_function_type)
        end_time = time.time()
        eval_time = end_time - start_time
        num_cont_pts = np.shape(control_points)[1]
        # print("control_points: " , control_points)
        print("computation time:" , end_time - start_time)
        spline_start_time = 0
        scale_factor = 1
        bspline = BsplineEvaluation(control_points, order, spline_start_time, scale_factor, False)
        number_data_points = 10000
        spline_data, time_data = bspline.get_spline_data(number_data_points)
        path_length = bspline.get_arc_length(number_data_points)
        waypoints = waypoint_data.get_waypoint_locations()
        ax[i,j].plot(spline_data[0,:], spline_data[1,:])
        ax[i,j].set_xlabel("x (m) \n evaluation time: " + str(np.round(eval_time,2)) + "\n num ctrl pts: " + str(num_cont_pts))
        ax[i,j].set_aspect('equal')
        plot2D_waypoints(waypoint_data, ax[i,j],arrow_scale = 2)
        if j == 0:
            alpha = 0.2
        else:
            alpha = 1
        plot_sfcs(sfc_list[1].get_sfc_list(), ax[i,j], alpha=alpha)

# curvature_data, time_data = bspline.get_spline_curvature_data(1000)
# acceleration_data, time_data = bspline.get_spline_derivative_data(1000,2)
# minvo_cps = bspline.get_minvo_control_points()
ax[0,0].set_title("No SFC's")
ax[0,1].set_title("With SFC Constraints")
ax[0,0].set_ylabel("Case 1 \n y (m)")
ax[1,0].set_ylabel("Case 2 \n y (m)")
ax[2,0].set_ylabel("Case 3 \n y (m)")
plt.show()

# plt.figure()
# plt.plot(time_data, curvature_data)
# plt.show()

# print("start curvature: " , curvature_data[0])
# print("end curvature: " , curvature_data[-1])


