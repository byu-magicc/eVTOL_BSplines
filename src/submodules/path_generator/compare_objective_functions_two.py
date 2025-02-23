
import numpy as np
import matplotlib.pyplot as plt
from bsplinegenerator.bsplines import BsplineEvaluation
# from path_generation.path_generator_compare import PathGenerator
import time
from path_generation.path_generator import PathGenerator
from path_generation.waypoint_data import Waypoint, WaypointData

## try different initial conditions
## move the max velocity around
## try direction constraint
## add function to get path length
max_curvature = 1.5
order = 3
waypoint_1 = Waypoint(location=np.array([[0],[0]]))
waypoint_2 = Waypoint(location=np.array([[6],[0]]))
dimension = np.shape(waypoint_1.location)[0]
# case = 1
# waypoint_1.velocity = np.array([[1],[0]])
# waypoint_2.velocity = np.array([[1],[0]])
# case = 2
# waypoint_1.velocity = np.array([[1],[0]])
# waypoint_2.velocity = np.array([[0],[1]])
# case = 3
# waypoint_1.velocity = np.array([[0],[1]])
# waypoint_2.velocity = np.array([[0],[-1]])
# case = 4
# waypoint_1.velocity = np.array([[-1],[0]])
# waypoint_2.velocity = np.array([[-1],[0]])
# case = 5
# waypoint_1.velocity = np.array([[0],[1]])
# waypoint_2.velocity = np.array([[0],[1]])
case = 6
waypoint_1.velocity = np.array([[-1],[0]])
waypoint_2.velocity = np.array([[1],[0]])

waypoint_data = WaypointData(start_waypoint=waypoint_1,end_waypoint=waypoint_2)

initial_control_points = None
# initial_control_points = np.array([[2,0,-2,-2,2.5,7,7,5,3],[0,0,0,4,4,4,0,0,0]]) # case 3 top
# initial_control_points = np.array([[2,0,-2,-2,2.5,7,7,5,3],[0,0,0,-4,-4,-4,0,0,0]]) # case 3 bottom
# initial_control_points = np.array([[5,-1,-2,3,2,7,6,0],[3,-2,3,4,-3,-3,1,-3]]) # case 3 middle
# initial_control_points = np.array([[0,0,0,2.5,2.5,2.5,5,5,5],[-3,0,3,3,0,-3,-3,0,3]]) # case 4
# initial_control_points = np.array([[2,0,-2,-2,2,2,3,5,7],[0,0,0,-2,-2,-2,0,0,0]]) # case 5

# curvature_method = "root_finding"
curvature_method = "analytical_roots"
# curvature_method = "control_point_derivatives"
# curvature_method = "constrain_max_acceleration_and_min_velocity"
# curvature_methods = ["root_finding", "analytical_roots", "control_point_hull", "control_point_derivatives_rotate", "constrain_max_acceleration_and_min_velocity"]
# curvature_methods = ["root_finding", "geometric","analytical_roots", "control_point_hull"]
# curvature_methods = ["geometric","analytical_roots", "control_point_hull", "constrain_max_acceleration_and_min_velocity"]
# curvature_methods = ["analytical_roots", "control_point_derivatives_rotate", "constrain_max_acceleration_and_min_velocity"]
# curvature_methods = ["analytical_roots","geometric","geometric","geometric"]
# colors = np.array(["r", "c"])
objective_methods = ["minimal_distance_path", "minimal_velocity_path",  "minimal_acceleration_path"]
line_styles = ['solid', 'dashed' , 'dotted']
hatches = [None, '/', '.']
num_methods = len(objective_methods)
color =  "c"
fig, ax = plt.subplots(3,num_methods)
max_x = 0
max_y = 0
min_x = 0
min_y = 0

for i in range(0,len(objective_methods)):
    path_gen = PathGenerator(dimension)
    start_time = time.time()
    control_points = path_gen.generate_path(waypoint_data=waypoint_data, max_curvature=max_curvature,
        max_incline=None, sfc_data=None, obstacles=None, objective_function_type=objective_methods[i])
    end_time = time.time()
    spline_start_time = 0
    scale_factor = 1
    bspline = BsplineEvaluation(control_points, order, spline_start_time, scale_factor, False)
    number_data_points = 1000000
    spline_data, time_data = bspline.get_spline_data(number_data_points)
    spline_at_knot_points, knot_points = bspline.get_spline_at_knot_points()
    curvature_data, time_data = bspline.get_spline_curvature_data(number_data_points)
    velocity_data, time_data = bspline.get_derivative_magnitude_data(number_data_points,1)
    velocity_std = np.std(velocity_data)
    # centripetal_acceleration_data, time_data = bspline.get_centripetal_acceleration_data(number_data_points)
    acceleration_magnitude_data, time_data = bspline.get_derivative_magnitude_data(number_data_points,2)
    ax[0,i].plot(spline_data[0,:],spline_data[1,:],label="path",color=color,linestyle=line_styles[i])
    ax[0,i].set_xlabel("x")
    # ax[0,i].set_ylabel("y")
    ax[1,i].plot(time_data,velocity_data,label="path",color=color,linestyle=line_styles[i])
    ax[1,i].set_xlabel("(sec)")
    max_x = np.max(spline_data[0,:])
    min_x = np.min(spline_data[0,:])
    max_y = np.max(spline_data[1,:])
    min_y = np.min(spline_data[1,:])
    width = max_x - min_x
    height = max_y - min_y
    ratio = abs((width)/(height))
    # ax[0,i].set_aspect(1/ratio)
    if case == 1:
        ax[0,i].set_ylim([-1, 1])
        ax[0,i].set_aspect((width)/(2))
    if objective_methods[i] == "minimal_distance_path":
        ax[0,i].set_title("Vel Cont Pt Mag",weight='bold')
    elif objective_methods[i] == "minimal_velocity_path":
        ax[0,i].set_title("Accel Cont Pt Mag",weight='bold')
    elif objective_methods[i] == "minimal_acceleration_path":
        ax[0,i].set_title("Jerk Cont Pt Mag",weight='bold')
    evaluation_time = end_time - start_time
    path_length = bspline.get_arc_length(number_data_points)
    acceleration_integral = np.sum(acceleration_magnitude_data)*time_data[1]
    curvature_integral = np.sum(curvature_data)*time_data[1]
    path_time = time_data[-1]
    # centripetal_acceleration_integral = np.sum(centripetal_acceleration_data)*time_data[1]/time_data[-1]
    # curvature_integral = np.sum(curvature_data)*time_data[1]
    curvature_extrema = np.max(curvature_data)
    ax[2,0].bar([i], [curvature_extrema] , color=color, hatch=hatches[i])
    ax[2,1].bar([i], [evaluation_time] , color=color, hatch=hatches[i])
    ax[2,2].bar([i], [path_length] , color=color, hatch=hatches[i])
    ax[2,i].tick_params(labelbottom = False, bottom = False)
    ax[2,0].set_ylabel('')
    ax[2,1].set_ylabel('(sec)')
    ax[2,2].set_ylabel('(m)')
    head_width = 0.3
    head_length = 0.3
    ax[0,i].arrow(waypoint_1.location.item(0), waypoint_1.location.item(1), waypoint_1.velocity.item(0)/2, 
                  waypoint_1.velocity.item(1)/2, head_width=head_width, head_length=head_length, color = 'k')
    ax[0,i].arrow(waypoint_2.location.item(0), waypoint_2.location.item(1), waypoint_2.velocity.item(0)/2, waypoint_2.velocity.item(1)/2,
                  head_width=head_width, head_length=head_length, color = 'k')
    ax[0,i].scatter(waypoint_1.location.item(0), waypoint_1.location.item(1),color='k', s=8)
    ax[0,i].scatter(waypoint_2.location.item(0), waypoint_2.location.item(1),color='k', s=8)


# for i in range(0,len(curvature_methods)):
#     ax[0,i].set_aspect(abs((max_x-min_x)/(max_y-min_y))*ratio)


ax[2,0].set_title("Max Curvature")
ax[2,1].set_title("Eval Time")
ax[2,2].set_title("Path Length")
ax[2,0].plot([-1,4],[max_curvature, max_curvature], color='k')
# ax[2,2].text(-0.5,max_curvature+0.15,"max curvature")
ax[0,0].set_ylabel("Optimized Paths \n y", weight='bold')
ax[1,0].set_ylabel("Velocity Mag \n (m/s)", weight='bold')
ax[2,0].set_ylabel("Metrics",weight='bold')
# ax[1,2].legend()

fig.supxlabel("Curvature Constraint Methods")
# fig.suptitle("Curvature Contrained Path Optimizations: Case " + str(case), weight='bold')
plt.tight_layout()
plt.subplots_adjust(wspace=0.5, hspace=0.5)
# plt.legend()
plt.show()

