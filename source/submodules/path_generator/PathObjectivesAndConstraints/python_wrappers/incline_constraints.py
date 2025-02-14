import ctypes 
import pathlib 
import os 
import numpy as np

script_dir = os.path.abspath(os.path.dirname(__file__))
libname_str = os.path.join(script_dir)
libname = pathlib.Path(libname_str)
lib = ctypes.CDLL(libname / "../build/src/libPathObjectivesAndConstraints.so")

class InclineConstraints(object):

    def __init__(self):
        self._order = 3
        ND_POINTER_DOUBLE = np.ctypeslib.ndpointer(dtype=np.float64, ndim=1,flags="C")
        nd_pointer_c_double = np.ctypeslib.ndpointer(dtype=ctypes.c_double)
        self._dimension = 3
        lib.InclinationConstraint_3.argtypes = [ctypes.c_void_p]
        lib.InclinationConstraint_3.restype = ctypes.c_void_p
        lib.get_interval_incline_constraints_3.argtypes = [ctypes.c_void_p, 
                ND_POINTER_DOUBLE, ctypes.c_int, ctypes.c_double, ctypes.c_double]
        lib.get_interval_incline_constraints_3.restype = nd_pointer_c_double
        lib.get_spline_incline_constraint_3.argtypes = [ctypes.c_void_p, 
                ND_POINTER_DOUBLE, ctypes.c_int, ctypes.c_double, ctypes.c_double]
        lib.get_spline_incline_constraint_3.restype = ctypes.c_double
        self.obj = lib.InclinationConstraint_3(0)
    
    def get_interval_incline_constraints(self, cont_pts, scale_factor, max_incline):
        num_cont_pts = np.shape(cont_pts)[1]
        cont_pts_array = cont_pts.flatten().astype('float64')
        num_intervals = num_cont_pts - self._order
        nd_pointer_c_double = np.ctypeslib.ndpointer(dtype=ctypes.c_double, shape=(num_intervals))
        lib.get_interval_incline_constraints_3.restype = nd_pointer_c_double
        constraints = lib.get_interval_incline_constraints_3(self.obj, cont_pts_array, num_cont_pts, scale_factor, max_incline)
        return constraints
    
    def get_spline_incline_constraint(self, cont_pts, scale_factor, max_incline):
        num_cont_pts = np.shape(cont_pts)[1]
        cont_pts_array = cont_pts.flatten().astype('float64')
        constraint = lib.get_spline_incline_constraint_3(self.obj, cont_pts_array, num_cont_pts, scale_factor, max_incline)
        return constraint/max_incline*10
    
# control_points = np.array([[ 5,  8, 6, 7, 3],
#                                [5, 11, 2, 7, 5],
#                                [8,  7, 5, 9, 2]])
# max_incline = 2
# scale_factor = 1.5
# inc_const = InclineConstraints()
# inc_constraints = inc_const.get_interval_incline_constraints(control_points, scale_factor, max_incline)
# print("inc_constraints: " , inc_constraints)

# inc_constraint = inc_const.get_spline_incline_constraint(control_points, scale_factor, max_incline)
# print("inc_constraint: " , inc_constraint)

# curvature_constraints = curve_const.get_interval_curvature_constraints(control_points, max_curvature)
# print("curvature_constraints: " , curvature_constraints)