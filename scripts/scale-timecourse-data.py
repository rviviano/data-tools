# Demean and enforce unit variance on timeseries data

import numpy as np

def scale_timecourse_data(V):
    """
    Demeans data assuming that each row is a timecourse and divides the row by 
    it's standard deviation to enforce unit variance.

    Input: Data where each row is a timecourse, therefore columns are timepoints

    Output: Scaled data in matrix format
    """    
    # Get the mean of each row
    V_mean = V.mean(axis=1)
    # Get the std of each row
    V_std = V.std(axis=1)
    # Change zeros to ones (avoid dividing zero columns by zero)
    V_std[V_std==0] = 1
    V_std = V_std.T[:, np.newaxis]
    # Feature scale the rows of V
    V = np.divide((V - V_mean[:, np.newaxis]), V_std)
    return V