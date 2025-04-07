"""
landing trajectory with minimum acceleration / jerk
        2/17/22 - RWB
        3/18/24 - RWB
"""
import numpy as np
from math import ceil
from scipy.interpolate import BSpline
from scipy.linalg import norm
from scipy.optimize import minimize
import matplotlib.pyplot as plt



class Planner:
    def __init__(self,
                 degree: int=2, 
                 num_segments: int=10,
                 t0: float=0.0,
                 tf: float=10.0,
                 ):
        self.degree = degree
        self.num_segment = num_segments
        self.t0 = t0
        self.tf = tf
        self.knots = self.set_knots_uniform(t0, tf, degree, num_segments)
        self.control_points = np.zeros((2,num_segments+degree+1)) 
        # spline class from scipy library
        self.spl = BSpline(t=self.knots, c=self.control_points, k=self.degree)

    def set_knots_uniform(self, 
                          t0: float, 
                          tf: float, 
                          degree: int, 
                          num_segments: int,
                          ):
        """ Create a uniformly spaced knot vector
        Args:
            t0 (float): start time for spline
            tf (float): stop time for spline
            degree (float): degree of spline
            num_segments (float): number of uniform time segments between t0 and tf
        Requires:
            num_segments > degree
        Returns:
            numpy array of knot points of length (num_segments + 2*degree)
        """
        if num_segments > degree:
            delta = (tf-t0)/num_segments
            knots = np.arange(start=t0-delta*degree, 
                              stop=tf+delta*(degree+1), 
                              step=delta)
            self.knots = knots
            self.t0 = t0
            self.tf = tf
            self.degree = degree
            self.num_segments = num_segments
            self.spl = BSpline(t=self.knots, c=self.ctrl_pts, k=self.degree+1)
        else:
            print("Error in Planner.knots_uniform():")
            print("The number of segments must be strictly greater than the degree.")

    def set_knots_clamped_uniform(t0, tf, degree, num_segments):
        """ Create a clamped uniformly spaced knot vector

        Args:
            t0 (float): start time for spline
            tf (float): stop time for spline
            degree (float): degree of spline
            num_segments (float): number of uniform time segments between t0 and tf

        Requires:
            num_segments > degree

        Returns:
            numpy array of knot points of length (num_segments + 2*degree)
        """
        if num_segments > degree:
            delta = (tf-t0)/num_segments
            knots = np.concatenate((
                t0 * np.ones(degree),
                np.concatenate((
                    np.arange(start=t0, stop=tf+delta, step=delta),
                    tf * np.ones(degree)),
                    axis=0)
                ), axis=0)
        else:
            print("Error in Planner.knots_clamped_uniform():")
            print("The number of segments must be strictly greater than the degree.")
            knots = []
        return knots
    
    def start_constraint(self, 
                         S: np.ndarray, % [start_position, start_velocity, ...]
                         )
        (n,ell) = S.shape()
        se
    
    def plot(self):
        ''' Plot the BSpline.  
            Plot both control points knot points 
        '''
        t0 = self.spl.t[0]  # first knot is t0
        tf = self.spl.t[-1]  # last knot is tf
        N = ceil((tf - t0)/0.01)  # number of points in time vector so spacing is 0.01
        t = np.linspace(t0, tf, N)  # time vector
        position = self.spl(t)
        # 3D trajectory plot
        fig = plt.figure(1)
        ax = fig.add_subplot(111, projection='3d')
        # plot spline (convert YX(-Z) -> NED)
        ax.plot(position[:, 1], position[:, 0], -position[:, 2],
                'b', label='spline')
        # plot control points (convert YX(-Z) -> NED)
        ax.plot(self.spl.c[:, 1], self.spl.c[:, 0], -self.spl.c[:, 2],
                '-o', label='control points')
        # plot knot points (convert YX(-Z) -> NED)
        pos_at_knots = self.spl(self.spl.t)
        ax.plot(pos_at_knots[:, 1], pos_at_knots[:, 0], -pos_at_knots[:, 2],
                color='k', marker='o', label='knot points')
        ax.legend()
        ax.set_xlabel('x', fontsize=16, rotation=0)
        ax.set_ylabel('y', fontsize=16, rotation=0)
        ax.set_zlabel('z', fontsize=16, rotation=0)
        ax.set_xlim3d([-10, 10])
        plt.show()

    def plot_vel_accel_jerk(self):
        ''' Plot the BSpline.  
            Plot both control points knot points 
        '''
        t0 = self.spl.t[0]  # first knot is t0
        tf = self.spl.t[-1]  # last knot is tf
        N = ceil((tf - t0)/0.01)  # number of points in time vector so spacing is 0.01
        t = np.linspace(t0, tf, N)  # time vector
        position = self.spl(t)
        vel, accel, jerk = splineDerivatives(self.spl, t)
        # 3D trajectory plot
        fig = plt.figure(1)
        ax = fig.add_subplot(111, projection='3d')
        # plot control points (convert YX(-Z) -> NED)
        ax.plot(self.spl.c[:, 1], self.spl.c[:, 0], -self.spl.c[:, 2],
                '-o', label='control points')
        # plot spline (convert YX(-Z) -> NED)
        ax.plot(position[:, 1], position[:, 0], -position[:, 2],
                'b', label='spline')
        ax.legend()
        ax.set_xlabel('x', fontsize=16, rotation=0)
        ax.set_ylabel('y', fontsize=16, rotation=0)
        ax.set_zlabel('z', fontsize=16, rotation=0)
        ax.set_xlim3d([-10, 10])
        # velocity plot
        fig = plt.figure(2)
        ax = fig.add_subplot(311)
        ax.plot(t, vel)
        ax.set_xlabel('time', fontsize=14, rotation=0)
        ax.set_ylabel('velocity', fontsize=14)
        # acceleration plot
        ax = fig.add_subplot(312)
        ax.plot(t, accel)
        ax.set_xlabel('time', fontsize=14, rotation=0)
        ax.set_ylabel('acceleration', fontsize=14)
        # jerk plot
        ax = fig.add_subplot(313)
        ax.plot(t, jerk, label='jerk')
        ax.set_xlabel('time', fontsize=14, rotation=0)
        ax.set_ylabel('jerk', fontsize=14)
        plt.show()



def controlPointsPosVel(p0, v0, pf, vf, knots, order, num_segments):
    N = num_segments + order
    ctrl_pts = np.zeros((N, 3))
    # set first two control points to specify initial position and velocity
    ctrl_pts[0, :] = p0.T[0,:]
    ctrl_pts[1,:] = p0.T + ( (knots[order + 1] - knots[1]) / order ) * v0.T
    # intermediate control points in straight line between p0 and pf
    for i in range(2, N-2):
        alpha = (i-1) / (num_segments+1)
        ctrl_pts[i,:] = (1-alpha) * p0.T + alpha * pf.T
    # set last two control points to specify initial position and velocity
    ctrl_pts[N-2,:] = pf.T - (knots[order + N-1]-knots[N-1])/order * vf.T
    ctrl_pts[N-1, :] = pf.T[0,:]
    return ctrl_pts

def landingTrajectory(p0, v0, pf, vf, t0, tf, num_segments=1):
    order = 3  # order of the spline
    knots = uniformKnots(t0, tf, order, num_segments)
    N = num_segments + order
    ctrl_pts = controlPointsPosVel(p0, v0, pf, vf, knots, order, num_segments)
    spl = BSpline(t=knots, c=ctrl_pts, k=order)
    x0 = np.reshape(ctrl_pts[2:N-2,:], (1,3*(N-4)))
    res = minimize(cost, x0, args=spl, method='nelder-mead', )
    spl.c[2:-2, :] = np.reshape(res.x, (N-4, 3))
    return spl


def cost(x, spl):
    t0 = spl.t[0]  # first knot is t0
    tf = spl.t[-1]  # last knot is tf
    N = ceil((tf - t0)/0.1)  # number of points in time vector so spacing is 0.1
    time = np.linspace(t0, tf, N)  # time vector
    m, n = np.shape(spl.c)
    num_intermediate_cpts = m-4
    spl.c[2:-2,:] = np.reshape(x, (num_intermediate_cpts, 3))
    vel, accel, jerk= splineDerivatives(spl, time)
    J = 1.0 * max(accel) + 1.0 * max(jerk)
    return J


def splineMagnitude(spl, time):
    mag = []
    for i in range(0, len(time)):
        mag.append(norm(spl(time[i])))
    return mag


def splineDerivatives(spl, time):
    spl_dot = spl.derivative()
    spl_ddot = spl_dot.derivative()
    spl_dddot = spl_ddot.derivative()
    vel = splineMagnitude(spl_dot, time)
    accel = splineMagnitude(spl_ddot, time)
    jerk = splineMagnitude(spl_dddot, time)
    return vel, accel, jerk




if __name__ == "__main__":
    plnr = Planner()
    #plnr.plot()
    #Planner.set_knots_clamped_uniform(t0=0, tf=1, degree=3, num_segments=10)
    
    # initial and final positions and velocities (NED)
    p0 = np.array([[100, 0, -10]]).T
    v0 = np.array([[-10, 0, 0]]).T
    pf = np.array([[0, 0, 0]]).T
    vf = np.array([[0, 0, 1]]).T
    # initial and final time
    t0 = 0
    tf = 5
    spl = landingTrajectory(p0, v0, pf, vf, t0, tf, num_segments=8)
    plotSpline(spl)

