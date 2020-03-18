# Stats utilities. Just a collection of convenience functions for basic
# descriptive stats or data cleaning.

# Raymond Viviano
# rayviviano@gmail.com
# January 2017

from __future__ import print_function
import os, sys, subprocess, traceback
import numpy as np
import pandas as pd
import scipy.stats
from scipy.stats import norm
from scipy.stats import gaussian_kde
import scipy.spatial.distance
import statsmodels.api as sm
from statsmodels.sandbox.regression.predstd import wls_prediction_std
from decimal import Decimal
import pprint
import time


### Function Definitions ###
def skewness(data):
    """Calculates skew by the adjusted Fisher-Pearson coeff for skew"""
    mu = np.nanmean(data)
    sigma = np.nanstd(data, ddof=1)
    n = float(np.count_nonzero(~np.isnan(data)))
    scaled_data = (data-mu)/sigma
    skew = (n/((n-1)*(n-2))) * np.nansum(scaled_data**3)
    return skew


def univariate_outlier_check(df, stdev_cutoff=3, skips=[]):
    """ Checks each variable in dataframe for univariate outliers 
        Skip columns by providing the column names in the 'skips' list.

        Flags datapoints 3 standard deviations away from the mean by default
    """
    colnames = df.columns.values
    df_idx = df.index.get_values()
    for col in colnames:
        if col not in skips:
            try:
                outlier_count = 0
                print('-'*80)
                print("Univariate outlier check for", col)
                np_col = df[col].as_matrix()
                z_col = (np_col - np.nanmean(np_col)) / np.nanstd(np_col)
                for case, z_score in zip(df_idx, z_col):
                    if (z_score > stdev_cutoff) or (z_score < -1*stdev_cutoff):
                        print(case, "z-score:", z_score)
                        outlier_count += 1
                if outlier_count == 0:
                    print("No univariate outliers exist for " + col + "...")
                else:
                    print("Total oulier count for " + col + ": ", + outlier_count)
            except:
                traceback.print_exc(file=sys.stdout)


def multivariate_outlier_check(df):
    """ Calculate Mahalanobis distance between one case and the centroid
        of all other cases in the dataset. Takes into count all relevant 
        variables. Make sure that the dataframe you pass only includes 
        variables that make sense as interval or ratio-scale numbers.
        
        TODO: Implement automatic checking of the Mahalanobis distance 
        off of a chi-squared distribution automatically determine outliers"""
    # Get dataframe row indices
    df_idx = df.index.get_values()

    # Convert dataframe to array
    df_mat  = df.as_matrix()

    # Calculate data centroid
    df_cent = np.mean(df_mat, axis=0)

    # Calculate data precision (inverse covariance) matrix
    df_prec = np.linalg.pinv(np.cov(df_mat.T))

    # List to hold all mahalanobis distances for future checking
    mahal_list = []

    print("Calculating Mahalanobis distances for multivariate outlier checking")
    for case, idx in zip(df_idx, [i for i in range(len(df_idx))]):
        case_vec = df.loc[case].as_matrix()
        mahal = scipy.spatial.distance.mahalanobis(case_vec, df_cent, df_prec)
        mahal_list.append(mahal)

    # Add mahalanobis distances to dataframe
    df['Mahalanobis'] = pd.Series(mahal_list, index=df.index)

    return df


def check_normality_assumptions(df):
    """ Check if a variable is approx normally disributed.
        Kinda lazy, just tries to print a bunch of basic descriptive stats
        and passes by expections"""

    print('-'*80)
    print('Normality assumption and data accuracy checks (skipping nans)')
    print('Note:Kurtosis is Fisher\'s (excess) kurtosis. 0 = Normal, not 3.')
    for i in df.columns.values:
        print('-'*80 + '\n' + i)
        try:   
            print("Min:",   np.nanmin(df[i]))
            print("Max:",   np.nanmax(df[i]))
            print("Mean: ", np.nanmean(df[i]))
            print("StDev: ", np.nanstd(df[i], ddof=1))
            print("Skewness: ", skewness(df[i].values))
            print("Kurtosis: ", df[i].kurtosis(skipna=True))        
        except:
            print('Error with ' + i)
            traceback.print_exc(file=sys.stdout)
            pass