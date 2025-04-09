"""
BSpline object
        3/26/25 - RWB
        4/1/25 - RWB
"""
import numpy as np
from math import floor, ceil
#from scipy.linalg import norm
import matplotlib.pyplot as plt

class BSpline:
    def __init__(self,
                 degree: int=2, # d
                 num_segments: int=10, # M
                 t0: float=0.0,
                 tf: float=10.0,
                 control_points: np.ndarray=np.random.random((3,10+2+1)),
                 ):
        self.degree = degree
        self.num_segments = num_segments
        self.t0 = t0
        self.tf = tf
        self.sigma_1 = self.num_segments / (self.tf - self.t0)
        self.sigma_0 = self.num_segments * self.t0 / (self.tf - self.t0)
        self.knots = self._create_uniform_knots()
        self.control_points = control_points
        self.n = control_points.shape[0]
        self.N_knot_res = 100  # resolution of the internal spline representation (1/N_knot_res)
        self.basis_dict = self._create_basis_dictionary()

    def _create_uniform_knots(self):
        ''' 
        Create a uniformly spaced knot vector
        Args:
            t0 (float): start time for spline
            tf (float): stop time for spline
            degree (float): degree of spline
            num_segments (float): number of uniform time segments between t0 and tf
        Requires:
            num_segments > degree
        Returns:
            numpy array of knot points of length (num_segments + 2*degree)
        '''
        if self.num_segments > self.degree:
            delta = (self.tf-self.t0)/self.num_segments
            knots = np.arange(start=self.t0-delta*self.degree, 
                              stop=self.tf+delta*(self.degree+1), 
                              step=delta)
        else:
            print('number of segments must exceed the degree')
        return knots

    def _create_basis_dictionary(self):
        '''
        returns a dictionary of numpy arrays where basis_dict['d','l'] is the basis vector for degree d, 
        l-th derivative
        '''
        M1 = np.array([[-1., 1.],
                       [1., 0.]])
        M2 = np.array([[1., -2., 1.],
                       [-2., 2., 1.],
                       [1., 0., 0.]])/2.
        M3 = np.array([[-2., 6., -6., 2.], 
                      [6., -12., 0., 8.],
                      [-6., 6., 6., 2.], 
                      [2., 0., 0., 0.]])/12.
        M4 = np.array([[1., -4., 6., -4., 1.], 
                       [-4., 12., -6., -12., 11.],
                       [6., -12., -6., 12., 11.], 
                       [-4., 4., 6., 4., 1.],
                       [1., 0., 0., 0., 0.]])/24.
        M5 = np.array([[-1., 5., -10., 10., -5., 1.],
                       [5., -20., 20., 20., -50., 26.],
                       [-10., 30., 0., -60., 0., 66.],
                       [10., -20., -20., 20., 50., 26.],
                       [-5., 5., 10., 10., 5., 1.],
                       [1., 0., 0., 0., 0., 0.]])/120.
        D0 = np.array([[-1.], 
                       [1.]])
        D1 = np.array([[-1., 0.], 
                       [1., -1.], 
                       [0., 1.]])   
        D2 = np.array([[-1., 0., 0.], 
                       [1., -1., 0.],
                       [0., 1., -1.],
                       [0., 0., 1.]])
        D3 = np.array([[-1., 0., 0., 0.], 
                       [1., -1., 0., 0.],
                       [0., 1., -1., 0.], 
                       [0., 0., 1., -1.],
                       [0., 0., 0., 1.]])
        D4 = np.array([[-1., 0., 0., 0., 0.], 
                       [1., -1., 0., 0., 0.],
                       [0., 1., -1., 0., 0.], 
                       [0., 0., 1., -1., 0.],
                       [0., 0., 0., 1., -1.],
                       [0., 0., 0., 0., 1.]])
        # create a dictionary of arrays
        eps = 1/self.N_knot_res
        t = np.arange(0., 1+eps, eps).reshape(1,self.N_knot_res+1)
        B = np.ones((1, t.shape[1]))
        basis_dict = {'d=0,l=0': B}
        for d in range(1, self.degree+1):
            tmp = B[0, :] * t
            B = np.concatenate((tmp, B), axis=0)
            if d==1:
                basis_dict['d=1,l=0'] = M1 @ B
                basis_dict['d=1,l=1'] = D0 @ basis_dict['d=0,l=0']
            elif d==2:
                basis_dict['d=2,l=0'] = M2 @ B
                basis_dict['d=2,l=1'] = D1 @ basis_dict['d=1,l=0']
                basis_dict['d=2,l=2'] = D1 @ D0 @ basis_dict['d=0,l=0']
            elif d==3:
                basis_dict['d=3,l=0'] = M3 @ B
                basis_dict['d=3,l=1'] = D2 @ basis_dict['d=2,l=0']
                basis_dict['d=3,l=2'] = D2 @ D1 @ basis_dict['d=1,l=0']
                basis_dict['d=3,l=3'] = D2 @ D1 @ D0 @ basis_dict['d=0,l=0']
            elif d==4:
                basis_dict['d=4,l=0'] = M4 @ B
                basis_dict['d=4,l=1'] = D3 @ basis_dict['d=3,l=0']
                basis_dict['d=4,l=2'] = D3 @ D2 @ basis_dict['d=2,l=0']
                basis_dict['d=4,l=3'] = D3 @ D2 @ D1 @ basis_dict['d=1,l=0']
                basis_dict['d=4,l=4'] = D3 @ D2 @ D1 @ D0 @ basis_dict['d=0,l=0']
            elif d==5:
                basis_dict['d=5,l=0'] = M5 @ B
                basis_dict['d=5,l=1'] = D4 @ basis_dict['d=4,l=0']
                basis_dict['d=5,l=2'] = D4 @ D3 @ basis_dict['d=3,l=0']
                basis_dict['d=5,l=3'] = D4 @ D3 @ D2 @ basis_dict['d=2,l=0']
                basis_dict['d=5,l=4'] = D4 @ D3 @ D2 @ D1 @ basis_dict['d=1,l=0']
                basis_dict['d=5,l=5'] = D4 @ D3 @ D2 @ D1 @ D0 @ basis_dict['d=0,l=0']
            else:
                print('Not implemented for d>5')
        return basis_dict
    
    def eval(self, 
             t, # time spline is being evaluated
             ell: int=0, # derivative being evaluated (ell=0 is original trajectory)
             ):
        '''
            Evaluate the spline at time t
                t can be a scalar or (N,) vector
        '''
        if type(t) is float:
            if (t<self.t0) or (t>self.tf):
                print('Error in BSpline.eval: t is not within limits')
            else:
                sigma = self.sigma_1 * t - self.sigma_0  # convert to 0 <= sigma <= M
                #m = int(np.floor(sigma))  # find interval [m, m+1] \subset [0, M]
                m = floor(sigma)  # find interval [m, m+1] \subset [0, M]
                ctrl = self.control_points[:, m:m+self.degree+1] # active control points c_{m:m+d}
                sigma = (self.sigma_1 * t - self.sigma_0) % 1.0  # convert to 0 <= sigma <= 1
                idx = round(sigma * self.N_knot_res)  # find index into basis array
                return ctrl @ self.basis_dict[f"d={self.degree},l={ell}"][:, idx:idx+1]
        else:
            if (t.item(0)<self.t0) or (t.item(-1)>self.tf):
                print('Error in BSpline.eval: t is not within limits')
            else:
                #---NOTE: Need to vectorize to make this more efficient
                traj = np.zeros((self.control_points.shape[0], t.shape[0]))
                for j in range(0, t.shape[0]-1):
                    sigma = self.sigma_1 * t[j] - self.sigma_0  # convert to 0 <= sigma <= M
                    m = int(np.floor(sigma))  # find interval [m, m+1] \subset [0, M]
                    ctrl = self.control_points[:, m:m+self.degree+1] # active control points c_{m:m+d}
                    sigma_ = sigma % 1.0  # convert to 0 <= sigma_ <= 1
                    idx = floor(sigma_ * self.N_knot_res)  # find index into basis array
                    traj[:,j:j+1] = ctrl @ self.basis_dict[f"d={self.degree},l={ell}"][:, idx:idx+1]
                return traj
    
    # def differentiate(self, 
    #                   ell: int=1):
    #     '''
    #         Returns a spline object that is the ell-th derivative of the current spline
    #     '''
    #     if ell > self.degree:
    #         print('cannot take more derivatives than the degree')
    #     else:
    #         C = self.control_points
    #         for l in range(1, ell+1):
    #             C[:,0:self.degree+self.num_segments-l+1] = C[:, 1:self.degree+self.num_segments-l+2] - C[:,0:self.degree+self.num_segments-l+1]
    #         diff = BSpline(
    #             degree=self.degree-ell,
    #             num_segments=self.num_segments,
    #             t0=self.t0,
    #             tf=self.tf,
    #             control_points=C[:,0:self.degree+self.num_segments-ell+1],
    #             )
    #     return diff
        
    def plot(self,
             fig_number: int=1,
             ):
        ''' 
        Plot the 3D BSpline.  
        Plot both control points knot points 
        '''
        N = ceil((self.tf - self.t0)/0.01)  # number of points in time vector so spacing is 0.01
        t = np.linspace(self.t0, self.tf, N)  # time vector
        position = self.eval(t)
        # 3D trajectory plot
        fig = plt.figure(fig_number, figsize=(6, int(f'{self.degree+6}')))
        ax = fig.add_subplot(2, 1, 1, projection='3d')
        # plot spline (convert YX(-Z) -> NED)
        ax.plot(position[1, :], position[0, :], -position[2, :],
                'b', label='spline')
        # plot control points (convert YX(-Z) -> NED)
        ax.plot(self.control_points[1, :], self.control_points[0, :], -self.control_points[2, :],
                color='g', marker='o', label='control points')
        # plot knot points (convert YX(-Z) -> NED)
        pos_at_knots = self.eval(self.knots[self.degree:self.degree+self.num_segments+1])
        ax.plot(pos_at_knots[1, :], pos_at_knots[0, :], -pos_at_knots[2, :],
                color='k', marker='o', linestyle='', label='knot points')
        ax.legend()
        ax.set_xlabel('x', fontsize=16, rotation=0)
        ax.set_ylabel('y', fontsize=16, rotation=0)
        ax.set_zlabel('z', fontsize=16, rotation=0)
        #---plot derivatives---
        for ell in range(1, self.degree+1):
            ax = fig.add_subplot(int(f'{self.degree+3}'), 1, int(f'{ell+3}'))
            tmp = self.eval(t, ell=ell)
            f = np.linalg.norm(tmp, axis=0)
            ax.plot(t, f)
            ax.set_xlabel('time', fontsize=10, rotation=0)
            ax.set_ylabel(f'ell={ell}', fontsize=10)
        plt.show()


if __name__ == "__main__":
    spl = BSpline(
        degree=3,  # eval bug for degree=4, plotting not quite right for degree=3
        num_segments=10, 
        t0=0.0,
        tf=10.0,
        control_points = np.array([[-2, -2, -2], 
                                   [-1, -1, -1],
                                   [0, 0, 0],
                                   [1, 1, 1],
                                   [2, 2, 2],
                                   [3, 3, 3],
                                   [4, 4, 4],
                                   [5, 5, 8],
                                   [6, 6, 6],
                                   [7, 7, 7],
                                   [8, 8, 8],
                                   [9, 9, 9],
                                   [10, 10, 10],
                                   ]).T
        )
    spl.plot(fig_number=1)


