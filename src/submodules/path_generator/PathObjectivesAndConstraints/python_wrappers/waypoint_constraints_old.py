import ctypes 
import pathlib 
import os 
import numpy as np

script_dir = os.path.abspath(os.path.dirname(__file__))
libname_str = os.path.join(script_dir)
libname = pathlib.Path(libname_str)
lib = ctypes.CDLL(libname / "../build/src/libPathObjectivesAndConstraints.so")

class WaypointConstraintsOld(object):

    def __init__(self, dimension):
        ND_POINTER_DOUBLE = np.ctypeslib.ndpointer(dtype=np.float64, ndim=1,flags="C")
        ND_POINTER_BOOL = np.ctypeslib.ndpointer(dtype=bool, ndim=1,flags="C")
        ND_POINTER_C_DOUBLE = np.ctypeslib.ndpointer(dtype=ctypes.c_double, shape=(dimension,2))
        ND_POINTER_C_DOUBLE_2 = np.ctypeslib.ndpointer(dtype=ctypes.c_double, shape=(2))
        self._dimension = dimension
        if dimension == 2:
            lib.WaypointConstraintsOld_2.argtypes = [ctypes.c_void_p]
            lib.WaypointConstraintsOld_2.restype = ctypes.c_void_p
            lib.velocity_at_waypoints_constraints_2.argtypes = [ctypes.c_void_p, 
                ND_POINTER_DOUBLE, ctypes.c_int, ctypes.c_double, ND_POINTER_DOUBLE, ND_POINTER_BOOL]
            lib.velocity_at_waypoints_constraints_2.restype = ND_POINTER_C_DOUBLE
            lib.acceleration_at_waypoints_constraints_2.argtypes = [ctypes.c_void_p, 
                ND_POINTER_DOUBLE, ctypes.c_int, ctypes.c_double, ND_POINTER_DOUBLE, ND_POINTER_BOOL]
            lib.acceleration_at_waypoints_constraints_2.restype = ND_POINTER_C_DOUBLE
            self.obj = lib.WaypointConstraintsOld_2(0)
        else: # value == 3
            lib.WaypointConstraintsOld_3.argtypes = [ctypes.c_void_p]
            lib.WaypointConstraintsOld_3.restype = ctypes.c_void_p
            lib.velocity_at_waypoints_constraints_3.argtypes = [ctypes.c_void_p, 
                ND_POINTER_DOUBLE, ctypes.c_int, ctypes.c_double, ND_POINTER_DOUBLE, ND_POINTER_BOOL]
            lib.velocity_at_waypoints_constraints_3.restype = ND_POINTER_C_DOUBLE
            lib.acceleration_at_waypoints_constraints_3.argtypes = [ctypes.c_void_p, 
                ND_POINTER_DOUBLE, ctypes.c_int, ctypes.c_double, ND_POINTER_DOUBLE, ND_POINTER_BOOL]
            lib.acceleration_at_waypoints_constraints_3.restype = ND_POINTER_C_DOUBLE
            self.obj = lib.WaypointConstraintsOld_3(0)

    def velocity_at_waypoints_constraints(self, cont_pts, scale_factor, desired_velocities, switches):
        num_cont_pts = np.shape(cont_pts)[1]
        cont_pts_array = cont_pts.flatten().astype('float64')
        des_vel_array = desired_velocities.flatten().astype('float64')
        switches_array = switches.flatten().astype('bool')
        if self._dimension == 2:
            objective = lib.velocity_at_waypoints_constraints_2(self.obj, cont_pts_array, num_cont_pts, scale_factor, des_vel_array, switches_array)
        else: # value = 3
            objective = lib.velocity_at_waypoints_constraints_3(self.obj, cont_pts_array, num_cont_pts, scale_factor, des_vel_array, switches_array)
        return objective
    
    def acceleration_at_waypoints_constraints(self, cont_pts, scale_factor, desired_accelerations, switches):
        num_cont_pts = np.shape(cont_pts)[1]
        cont_pts_array = cont_pts.flatten().astype('float64')
        des_accel_array = desired_accelerations.flatten().astype('float64')
        switches_array = switches.flatten().astype('bool')
        if self._dimension == 2:
            objective = lib.acceleration_at_waypoints_constraints_2(self.obj, cont_pts_array, num_cont_pts, scale_factor, des_accel_array, switches_array)
        else: # value = 3
            objective = lib.acceleration_at_waypoints_constraints_3(self.obj, cont_pts_array, num_cont_pts, scale_factor, des_accel_array, switches_array)
        return objective

    
# control_points = np.array([[3, 0, 9, 5, 7, 8, 0, 1],
#                             [1, 6, 2, 2, 1, 4, 1, 4]])
# desired_velocities = np.array([[7, 3],
#                                [1, 4]])
# switches = np.array([True, True])
# ### desired_velocities = desired_velocities / np.linalg.norm(desired_velocities,2,0)
# scale_factor = 1
# const_func = WaypointConstraintsOld(2)
# vel_at_waypoints_constraints = const_func.velocity_at_waypoints_constraints(control_points, scale_factor, desired_velocities, switches)
# print("vel_at_waypoints_constraints: " , vel_at_waypoints_constraints)

# desired_accelerations = np.array([[4, 4],
#                                [2, 7]])
# accel_at_waypoints_constraints = const_func.acceleration_at_waypoints_constraints(control_points, scale_factor, desired_accelerations, switches)
# print("accel_at_waypoints_constraints: " , accel_at_waypoints_constraints)


# control_points = np.array([[3, 0, 9, 5, 7, 8, 0, 1],
#                             [1, 6, 2, 2, 1, 4, 1, 4],
#                             [9, 4, 9, 6, 7, 12, 0, 2]])
# desired_velocities = np.array([[3, 3],
#                                [1, 1],
#                                [2, 7]])
# switches = np.array([True, True])
# ### desired_velocities = desired_velocities / np.linalg.norm(desired_velocities,2,0)
# scale_factor = 1
# const_func = WaypointConstraintsOld(3)
# vel_at_waypoints_constraints = const_func.velocity_at_waypoints_constraints(control_points, scale_factor, desired_velocities, switches)
# print("vel_at_waypoints_constraints: " , vel_at_waypoints_constraints)

# desired_accelerations = np.array([[4, 4],
#                                [2, 5],
#                                [1, 1]])
# accel_at_waypoints_constraints = const_func.acceleration_at_waypoints_constraints(control_points, scale_factor, desired_accelerations, switches)
# print("accel_at_waypoints_constraints: " , accel_at_waypoints_constraints)