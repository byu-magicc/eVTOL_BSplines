#this file is a helper file for my own benefit to understand how to work with the SFC constraints with David.
#At this point, I am frankly quite lost.


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


order = 3
dimension = 2
num_intervals_free_space = 5
# objective_function_type="minimal_acceleration_path"
objective_function_type = "minimal_distance_path"
max_curvature=0.5

point_1 = np.array([[0],[0]])
point_2 = np.array([[5],[0]])
point_3 = np.array([[5],[10]])
point_4 = np.array([[-1],[10]])
point_sequence = np.concatenate((point_1,point_2,point_3,point_4),axis=1)
R1, T1, min_len_1 = get2DRotationAndTranslationFromPoints(point_1, point_2)
R2, T2, min_len_2 = get2DRotationAndTranslationFromPoints(point_2, point_3)
R3, T3, min_len_3 = get2DRotationAndTranslationFromPoints(point_3, point_4)

#creates the corresponding safe flight corridors
sfc_1 = SFC(np.array([[min_len_1+4],[3]]), T1, R1)
sfc_2 = SFC(np.array([[min_len_2+3],[4]]), T2, R2)
sfc_3 = SFC(np.array([[min_len_3+4],[3]]), T3, R3)
sfcs = [sfc_1, sfc_2, sfc_3]
sfc_data = SFC_Data(sfcs,point_sequence)

# waypoints case 1
waypoint_a = Waypoint(location=point_1, velocity = (point_2-point_1)/np.linalg.norm((point_2-point_1)))
waypoint_b = Waypoint(location=point_4, velocity = (point_4-point_3)/np.linalg.norm((point_4-point_3)))
waypoint_data = WaypointData(start_waypoint=waypoint_a, end_waypoint=waypoint_b)



#creates the path generator
path_gen = PathGenerator(dimension=dimension,
                         num_intervals_free_space=num_intervals_free_space)


#calls the function to get the sfc constraints for this particular thing
sfc_constriants = path_gen._PathGenerator__create_safe_flight_corridor_constraint(sfc_data=sfc_data,
                                                                    num_cont_pts=7,
                                                                    num_intermediate_waypoints=0)




'''
start_time = time.time()
control_points, status =\
      path_gen.generate_path(waypoint_data=waypoint_data,
                             max_curvature=max_curvature,
                             max_incline=None,
                             sfc_data=sfc_data,
                             obstacles=None,
                             objective_function_type=objective_function_type)

end_time = time.time()
evaluation_time = end_time - start_time


num_cont_pts = np.shape(control_points)[1]
print("computation time: ", evaluation_time)

bspline = BsplineEvaluation(control_points=control_points,
                            order=order,
                            start_time=0.0,
                            scale_factor=1.0,
                            clamped=False)

num_data_points = 10000
spline_data, time_data = bspline.get_spline_data(num_data_points_per_interval=num_data_points)
path_length = bspline.get_arc_length(resolution=num_data_points)

waypoints = waypoint_data.get_waypoint_locations()

fix, ax = plt.subplots(1,1)
ax.plot(spline_data[0,:], spline_data[1,:])
ax.set_xlabel("x (m) \n evaluation time: " + str(np.round(evaluation_time,2)) + "\n num ctrl pts: " + str(num_cont_pts))
ax.set_aspect('equal')
plot2D_waypoints(waypoint_data, ax,arrow_scale = 2)
plot_sfcs(sfcs=sfcs, ax=ax, alpha=1.0)

plt.show()
#'''



potato = 0