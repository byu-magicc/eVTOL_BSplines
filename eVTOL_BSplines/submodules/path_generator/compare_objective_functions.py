
import numpy as np
import matplotlib.pyplot as plt
from bsplinegenerator.bsplines import BsplineEvaluation
from path_generation.path_generator import PathGenerator
from path_generation.waypoint_data import Waypoint, WaypointData
import time

## try different initial conditions
## move the max velocity around
## try direction constraint
## add function to get path length
max_curvature = 1.5
order = 3
waypoints = np.array([[0,6],[0,0]])
dimension = np.shape(waypoints)[0]
velocities_1 = np.array([[1,1],[0,0]]) # 1
velocities_2 = np.array([[1,0],[0,1]]) # 2
velocities_3 = np.array([[0,0],[1,-1]]) # 3
velocities_4 = np.array([[-1,-1],[0,0]]) # 4
velocities_5 = np.array([[0,0],[1,1]]) # 5
velocities_6 = np.array([[-1,1],[0,0]]) # 6
velocities_1 = velocities_1/np.linalg.norm(velocities_1,2,0) 
velocities_2 = velocities_2/np.linalg.norm(velocities_2,2,0) 
velocities_3 = velocities_3/np.linalg.norm(velocities_3,2,0)
velocities_4 = velocities_4/np.linalg.norm(velocities_4,2,0)
velocities_5 = velocities_5/np.linalg.norm(velocities_5,2,0)
velocities_6 = velocities_6/np.linalg.norm(velocities_6,2,0)
cases = (velocities_1, velocities_2, velocities_3, velocities_4, velocities_5, velocities_6)

path_length_averages = np.zeros(3)
eval_time_averages = np.zeros(3)
vel_sd_averages = np.zeros(3)


curvature_method = "analytical_roots"
objective_methods = ["minimal_distance_path", "minimal_velocity_path",  "minimal_acceleration_path"]
labels = ["velocity \n ctrl pt mag", "acceleration \n ctrl pt mag",  "jerk \n ctrl pt mag"]
line_styles = ['solid', 'dashed' , 'dotted']
hatches = [None, '/', '.']
color =  "c"


for i in range(len(objective_methods)):
    path_length_sum = 0
    eval_time_sum = 0
    vel_sd_sum = 0
    for j in range(len(cases)):
        path_gen = PathGenerator(dimension)
        velocities = cases[j]
        start_time = time.time()
        waypoint_1 = Waypoint(location=np.array([[0],[0]]))
        waypoint_2 = Waypoint(location=np.array([[6],[0]]))
        waypoint_1.velocity = velocities[:,0]
        waypoint_2.velocity = velocities[:,1]
        waypoint_data = WaypointData(start_waypoint=waypoint_1,end_waypoint=waypoint_2)
        start_time = time.time()
        control_points = path_gen.generate_path(waypoint_data=waypoint_data, max_curvature=max_curvature,
            max_incline=None, sfc_data=None, obstacles=None, objective_function_type=objective_methods[i])
        end_time = time.time()
        spline_start_time = 0
        scale_factor = 1
        bspline = BsplineEvaluation(control_points, order, spline_start_time, scale_factor, False)
        number_data_points = 1000000
        spline_data, time_data = bspline.get_spline_data(number_data_points)
        velocity_data, time_data = bspline.get_derivative_magnitude_data(number_data_points,1)
        velocity_std = np.std(velocity_data)
        curvature_data, time_data = bspline.get_spline_curvature_data(number_data_points)
        evaluation_time = end_time - start_time
        path_length = bspline.get_arc_length(number_data_points)
        curvature_extrema = np.max(curvature_data)
        path_length_sum += path_length
        eval_time_sum += evaluation_time
        vel_sd_sum += velocity_std
    path_length_averages[i] = path_length_sum/len(cases)
    eval_time_averages[i] = eval_time_sum/len(cases)
    vel_sd_averages[i] = vel_sd_sum/len(cases)
   

fig, ax = plt.subplots(1,3)

ax[0].set_title("Path Length", weight="bold")
ax[0].set_ylabel('length (m)')
ax[1].set_title("Evaluation Time", weight="bold")
ax[1].set_ylabel('time (sec)')
ax[2].set_title("Velocity Standard Deviation", weight="bold")
ax[2].set_ylabel('vel sd (m/s)')
for i in range(3):
    ax[0].bar([labels[i]], [path_length_averages[i]] , color=color, hatch=hatches[i])
    ax[1].bar([labels[i]], [eval_time_averages[i]] , color=color, hatch=hatches[i])
    ax[2].bar([labels[i]], [vel_sd_averages[i]] , color=color, hatch=hatches[i])



fig.supxlabel("Objective Functions")
# fig.suptitle("Curvature Contrained Path Optimizations: Case " + str(case), weight='bold')
plt.tight_layout()
plt.subplots_adjust(wspace=0.2, hspace=0.1)
# plt.legend()
plt.show()

