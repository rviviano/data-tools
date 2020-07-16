import numpy as np
import scipy.spatial.distance

def kmeans(data, centers, delta=.0001, maxiter=500, metric="cityblock"):
    """ 
    Inputs: 
        data    : observations x features
        centers : k x features: the initial centroids; e.g. 
                  a random sample of rows from X. These could be found 
                  using k++ initialization before calling kmeans.
        delta   : cutoff; iterate until the average distance of the 
                  observations to the centeroids is within delta of the 
                  previous average distance to the centroids; i.e, 
                  prevavgdist*(1-delta) <= current avgdist <= prevavgdist
        maxiter : maximum iterations of the algorithm if the
                  delta cutoff isn't reached
        metric  : cityblock = l1 norm or manhattan distance
        
    Outputs:
        centroids : k x features; the final centroids
        xlbl      : codebook; maps each feature in data to its nearest centroid, k 
    """        
    N, dim = data.shape
    k = centers.shape[0]         
    all_obs = np.arange(N)
    prevdist = 0
    for jth_iter in range(maxiter):
        D = scipy.spatial.distance.cdist(data, centers, metric=metric)
        xlbl = D.argmin(axis=1)  
        distances = D[all_obs, xlbl]
        avdist = distances.mean()
        
        if ((1-delta)*prevdist <= avdist <= prevdist) or (jth_iter == maxiter):
            break
            
        prevdist = avdist
        
        # Update centers
        for jc in range(k):
            c = np.where(xlbl == jc)[0]
            if len(c) > 0:
                centers[jc] = data[c].mean(axis=0)
    return centers, xlbl, jth_iter




