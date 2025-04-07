"""
generate uniform canonical basis 
messing with b-splines
        2/18/21 - RWB
        10/19/22 - RWB
"""
import numpy as np
from scipy.interpolate import BSpline
import matplotlib.pyplot as plt


def uniformSplineBasis(degree=0, plot=False):
        M = degree+1
        knots = list(range(-degree, M+degree+1))
        print(knots)
        controlPoints = np.zeros((M+degree,1))
        controlPoints[degree] = 1.0
        spl = BSpline(knots, controlPoints, degree)
        N = degree+1
        time = np.linspace(0, M, N*M+1)
        basisPoints = spl(time)
        if plot is True:
                plotFunction(time, basisPoints)
        return time, basisPoints 


def plotFunction(t, y):
        fig = plt.figure(0)
        plt.plot(t, y)


def plotSplineBasis(degree, M):
        fig = plt.figure(degree)
        t = np.linspace(0, M, 100)
        knots = uniformKnots(degree, M)
        fig.suptitle(f"degree={degree}, knots = {str(knots)}")
        ax = fig.subplots(M+degree)
        for i in range(0, M+degree):
                ctrl = [0] * (M+degree)
                ctrl[i] = 1
                pts = splineBasis(degree, knots, ctrl, t)
                ax[i].plot(t, pts)
                ax[i].set(ylabel=f"m={i}")

if __name__ == "__main__":

        plt.show()

        savestr = ''
        for degree in range(0, 8):
                foo = 'time' + str(degree) + ', basisPoints' + str(degree) + '= uniformSplineBasis(degree=degree)'
                exec(foo)
                savestr = savestr + 'time' + str(degree) + '=time' + str(degree) + ','
                savestr = savestr + 'basisPoints' + str(degree) + '=basisPoints' + str(degree)
                if degree<7:
                        savestr = savestr + ','
        foo = "np.savez('uniformSplineBasis'," + savestr + ")"
        exec(foo) 



