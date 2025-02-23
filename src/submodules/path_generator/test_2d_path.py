import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from bsplinegenerator.bsplines import BsplineEvaluation
# from path_generation.path_generator import PathGenerator
from path_generation.path_generator import PathGenerator
from path_generation.safe_flight_corridor import SFC, SFC_Data, plot_sfcs, get2DRotationAndTranslationFromPoints
from path_generation.obstacle import Obstacle, plot_2D_obstacles
from path_generation.waypoint_data import Waypoint, WaypointData
from path_generation.path_plotter import set_axes_equal
import time

point_1 = np.array([[0],[0]])
point_2 = np.array([[6],[0]])

waypoint_1 = Waypoint(location=point_1)
waypoint_2 = Waypoint(location=point_2)
waypoint_1.velocity = np.array([[1],[0]])
waypoint_2.velocity = np.array([[-1],[0]])
waypoint_data = WaypointData(start_waypoint=waypoint_1,end_waypoint=waypoint_2)


waypoints = waypoint_data.get_waypoint_locations()

max_curvature = 0.5
max_incline = None
order = 3

dimension = np.shape(point_1)[0]


path_gen = PathGenerator(dimension)
start_time = time.time()
control_points, status = path_gen.generate_path(waypoint_data=waypoint_data, max_curvature=max_curvature,
    max_incline=max_incline, sfc_data=None, obstacles=None)
end_time = time.time()
# print("control_points: " , control_points)
print("computation time:" , end_time - start_time)
spline_start_time = 0
scale_factor = 1
bspline = BsplineEvaluation(control_points, order, spline_start_time, scale_factor, False)
number_data_points = 10000
spline_data, time_data = bspline.get_spline_data(number_data_points)

spline_data, time_data = bspline.get_spline_data(1000)
curvature_data, time_data = bspline.get_spline_curvature_data(1000)
# acceleration_data, time_data = bspline.get_spline_derivative_data(1000,2)
minvo_cps = bspline.get_minvo_control_points()

# print("sfcs: " , sfcs)
ax = plt.axes()
ax.plot(spline_data[0,:], spline_data[1,:])
ax.scatter(waypoints[0,:],waypoints[1,:])
set_axes_equal(ax, dimension)
plt.show()

# plt.figure()
# plt.plot(time_data, curvature_data)
# plt.show()

print("start curvature: " , curvature_data[0])
print("end curvature: " , curvature_data[-1])


