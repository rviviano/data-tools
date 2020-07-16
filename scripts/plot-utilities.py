# A few plotting utilities for 2d & 3d scatterplots, and matrices

# Raymond Viviano
# July 28th, 2016

import numpy as np
import scipy.linalg
import statsmodels.api as sm
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

def best_fit_2nd_order_plane(data):
    """ Best fit plane (2nd order polynomial) for provided data points. Use to
        plot a plane in a 3-d plot """
    # Regular grid covering the domain of the data
    mn = np.min(data, axis=0)
    mx = np.max(data, axis=0)
    X,Y = np.meshgrid(np.linspace(mn[0]-2, mx[0]+2, 20), np.linspace(mn[1]-2, mx[1]+2, 20))
    XX = X.flatten()
    YY = Y.flatten()
    # Best-fit quadratic curve
    A = np.c_[np.ones(data.shape[0]), data[:,:2], np.prod(data[:,:2], axis=1), data[:,:2]**2]
    C,_,_,_ = scipy.linalg.lstsq(A, data[:,2])
    # Evaluate it on a grid
    Z = np.dot(np.c_[np.ones(XX.shape), XX, YY, XX*YY, XX**2, YY**2], C).reshape(X.shape)
    return X, Y, Z


def fit_line(x, y):
    """ Return slope, intercept of best fit line. """
    X = sm.add_constant(x)
    model = sm.OLS(y, X, missing='drop')
    fit = model.fit()
    return fit.params[1], fit.params[0]


def plot_matrix(m, ylabel=""):
    """ Simple way to plot a connectivity or var-covar matrix"""
    abs_max = abs(m).max()
    plt.imshow(m, cmap=plt.cm.RdBu_r, interpolation="nearest",
              vmin=-abs_max, vmax=abs_max)