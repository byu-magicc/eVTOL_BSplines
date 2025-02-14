import ctypes 
import pathlib 
import os 
import numpy as np

script_dir = os.path.abspath(os.path.dirname(__file__))
libname_str = os.path.join(script_dir)
libname = pathlib.Path(libname_str)
lib = ctypes.CDLL(libname / "../build/src/libPathObjectivesAndConstraints.so")

class WaypointConstraints(object):

    def __init__(self, dimension):
        ND_POINTER_DOUBLE = np.ctypeslib.ndpointer(dtype=np.float64, ndim=1,flags="C")
        ND_POINTER_C_DOUBLE = np.ctypeslib.ndpointer(dtype=ctypes.c_double, shape=(dimension))
        self._dimension = dimension
        if dimension == 2:
            lib.WaypointConstraints_2.argtypes = [ctypes.c_void_p]
            lib.WaypointConstraints_2.restype = ctypes.c_void_p
            lib.velocity_constraint_2.argtypes = [ctypes.c_void_p, 
                ND_POINTER_DOUBLE, ctypes.c_int, ND_POINTER_DOUBLE, ctypes.c_double, ctypes.c_bool]
            lib.velocity_constraint_2.restype = ND_POINTER_C_DOUBLE
            lib.acceleration_constraint_2.argtypes = [ctypes.c_void_p, 
                ND_POINTER_DOUBLE, ctypes.c_int, ND_POINTER_DOUBLE, ctypes.c_double, ctypes.c_bool]
            lib.acceleration_constraint_2.restype = ND_POINTER_C_DOUBLE
            self.obj = lib.WaypointConstraints_2(0)
        else: # value == 3
            lib.WaypointConstraints_3.argtypes = [ctypes.c_void_p]
            lib.WaypointConstraints_3.restype = ctypes.c_void_p
            lib.velocity_constraint_3.argtypes = [ctypes.c_void_p, 
                ND_POINTER_DOUBLE, ctypes.c_int, ND_POINTER_DOUBLE, ctypes.c_double, ctypes.c_bool]
            lib.velocity_constraint_3.restype = ND_POINTER_C_DOUBLE
            lib.acceleration_constraint_3.argtypes = [ctypes.c_void_p, 
                ND_POINTER_DOUBLE, ctypes.c_int, ND_POINTER_DOUBLE, ctypes.c_double, ctypes.c_bool]
            lib.acceleration_constraint_3.restype = ND_POINTER_C_DOUBLE
            self.obj = lib.WaypointConstraints_3(0)

    def velocity_constraint(self, cont_pts, desired_velocity,inverse_scale_factor, isStartVelocity):
        num_cont_pts = np.shape(cont_pts)[1]
        cont_pts_array = cont_pts.flatten().astype('float64')
        desired_velocity_array = desired_velocity.flatten().astype('float64')
        
        if self._dimension == 2:
            constraint = lib.velocity_constraint_2(self.obj, cont_pts_array, num_cont_pts, desired_velocity_array, inverse_scale_factor, isStartVelocity)
        else: # value = 3
            constraint = lib.velocity_constraint_3(self.obj, cont_pts_array, num_cont_pts, desired_velocity_array, inverse_scale_factor, isStartVelocity)
        return constraint
    
    def acceleration_constraint(self, cont_pts, desired_acceleration,inverse_scale_factor, isStartVelocity):
        num_cont_pts = np.shape(cont_pts)[1]
        cont_pts_array = cont_pts.flatten().astype('float64')
        desired_acceleration_array = desired_acceleration.flatten().astype('float64')
        if self._dimension == 2:
            constraint = lib.acceleration_constraint_2(self.obj, cont_pts_array, num_cont_pts, desired_acceleration_array, inverse_scale_factor, isStartVelocity)
        else: # value = 3
            constraint = lib.acceleration_constraint_3(self.obj, cont_pts_array, num_cont_pts, desired_acceleration_array, inverse_scale_factor, isStartVelocity)
        return constraint

    
# cont_pts = np.array([[3, 0, 9, 5, 7, 8, 0, 1],
#                             [1, 6, 2, 2, 1, 4, 1, 4]])
# desired_velocity = np.array([[7],
#                                [1]])
# is_start = True
# # ### desired_velocities = desired_velocities / np.linalg.norm(desired_velocities,2,0)
# inverse_scale_factor = 1
# const_func = WaypointConstraints(2)
# vel_at_waypoints_constraints = const_func.velocity_constraint(cont_pts, desired_velocity,inverse_scale_factor, is_start)
# print("vel_at_waypoints_constraints: " , vel_at_waypoints_constraints)

# cont_pts = np.array([[-3,  -4, -2, -.5, 1  ,   0,  2, 3.5],
#                     [.5, 3.5,  6, 5.5, 3.7,   2, -1,  2],
#                     [1, 3.2,  5,   0, 3.3, 1.5, -1, 2.5]])
# const_func = WaypointConstraints(3)
# inverse_scale_factor = 1/2
# is_start = True
# desired_acceleration = np.array([[0],
#                                [0],
#                                [0]])
# accel_at_waypoints_constraints = const_func.acceleration_constraint(cont_pts, desired_acceleration,inverse_scale_factor, is_start)
# print("accel_at_waypoints_constraints: " , accel_at_waypoints_constraints)