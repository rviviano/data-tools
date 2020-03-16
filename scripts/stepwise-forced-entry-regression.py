# Stepwise forced-entry regression

# Raymond Viviano
# rayviviano@gmail.com
# July 24th, 2017

from __future__ import print_function
import os, sys, subprocess, traceback
import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm


### Function Definitions ###
def stepwise_forced_entry(dv, ivs, order, df):
    """
    Performs forced entry stepwise regression.
    
    Arguments: 
        dv      : <str>                 Dependent variable column name
        
        ivs     : <list> <str>          List of independent var column names
        
        order   : <list> <int>          Order of iv entry into model, length 
                                        must equal ivs length.
                                        
        df      : <pandas dataframe>    Data for the regression 
    
    Returns: 
        prints all important info to standard output
    
    Example Function Call:
        c = stepwise_forced_entry(dv='HC_PCC_FC', ivs=['Age','Gender','SMC'],
                                  order=[0,0,1], df=alzheimers_data)
    
    Example Output:
        Step 0: 'HC_PCC_FC ~ Age + Gender'
        Step 1: 'HC_PCC_FC ~ Age + Gender + SMC'
        
        Step 0:
            R-sqr = ????, Adj-R2-sqr = ????,
            F-val = ????, p-val = ????, 
            N = ??, DF Resid = ??, DF Model = ??    
            -------------------------------------------------
                        coef        stderr      t       p
            Intercept   ????        ????        ????    ????
            Age         ????        ????        ????    ????
            Gender      ????        ????        ????    ????
            
        Step 1:
            R-sqr = ????, Adj-R-sqr = ????
            F-val = ????, p-val = ????, 
            N = ??, DF Resid = ??, DF Model = ??    
            -------------------------------------------------
                        coef        stderr      t       p
            Intercept   ????        ????        ????    ????
            Age         ????        ????        ????    ????
            Gender      ????        ????        ????    ????
            SMC         ????        ????        ????    ????
        
        Step 0 --> Step 1 Delta R-sqrd  =   ????    
        Step 0 --> Step 1 F             =   ????
        Step 0 --> Step 1 p-val         =   ????
        
    Notes:
        Uses listwise deletion

        Mod0 = Model for Step 0; Mod1 = Model for Step 1
        kMod = Model Degrees of Freedom = ivs in model
        N = Number of observations/data points
        
        Delta Rsqr F = ((Mod1R2-Mod0R2)/(kMod1-kMod0))/((1-Mod1R2)/(N-kMod1-1))
        dfn          = kMod1 - kMod0
        dfd          = N - kMod1 - 1
        p            = scipy.stats.f.sf(f,dfn,dfd)
        
    """
    # Separate each step with a line of hyphens for easier reading
    print('-'*80)

    # Sort lists together based on iv entry order 
    order, ivs = (list(x) for x in zip(*sorted(zip(order,ivs))))
        
    # Subset input dataframe to only use columns of interest
    all_vars = [x for x in ivs]
    all_vars.append(dv)
    df_sub = df.loc[:, all_vars]
    
    # Drop rows with missing data (Listwise Deletion)
    df_sub = df_sub.dropna()

    # Remove '-' and '.' from column names, dv, and ivs
    # Required for statsmodels formulas to run properly
    colnames = df_sub.columns.values
    new_colnames = [x.replace('-','_') for x in colnames]
    new_colnames = [x.replace('.','_') for x in new_colnames]
    df_sub.columns = new_colnames
    dv = dv.replace('-','_').replace('.','_')
    ivs = [x.replace('-','_').replace('.','_') for x in ivs]
    
    # Dictionary to hold results of each step
    step_results = {}
    
    # Stepwise forced entry regression
    for step in range(np.min(order), np.max(order)+1):
        # Devise statsmodels formula
        frm = dv + ' ~ '
        for i in range(len(order)):
            if order[i] <= step:
                frm += ivs[i] + ' + '
        frm = frm[:-2] # get rid of the extra plus from the above loop
        
        print('Step ' + str(step) + ': ' + frm)
        
        # Unstandardized regression
        model = sm.OLS.from_formula(formula=frm, data=df_sub)
        results = model.fit()
        step_results[step] = results
        nobs = len(df_sub[dv])
        
        print('R-sqr = ' + str(np.round(results.rsquared, 4)) + 
               ',\tAdj-R2-sqr = ' + str(np.round(results.rsquared_adj, 4)) + 
               '\nF-val = ' + str(np.round(results.fvalue, 4)) + ',\tp-val = '+ 
               str(np.round(results.f_pvalue, 4)) + '\nN = ' + str(nobs) + 
               ', DF Resid = ' + str(results.df_resid) + ', DF Model = ' + 
               str(results.df_model))
                       
        print('')

        df_results = pd.DataFrame({'coef': results.params, 
                                   'stdErr':results.bse, 
                                   't':results.tvalues, 
                                   'p':results.pvalues})
        # Reorder columns
        df_results = df_results[['coef', 'stdErr', 't', 'p']]

        
        print(df_results)  
        print('')

        # Change in R2 stats
        if step != np.min(order):
            N = nobs
            Mod0R2 = step_results[step -1].rsquared
            Mod1R2 = step_results[step].rsquared
            kMod0 = step_results[step -1].df_model
            kMod1 = step_results[step].df_model
            dfn = kMod1 - kMod0
            dfd = N - kMod1 - 1
            delta_rsqr = Mod1R2 - Mod0R2
            delta_rsqr_f = ((Mod1R2-Mod0R2)/(kMod1-kMod0))/((1-Mod1R2)/(N-kMod1-1))
            delta_rsqr_p = stats.f.sf(delta_rsqr_f,dfn,dfd)
            print('\nStep ', step-1, ' --> Step', step, 'Delta R-Sqrd     = ', delta_rsqr)    
            print('Step ', step-1, ' --> Step', step, 'Delta R-Sqrd F-Val = ', delta_rsqr_f)
            print('Step ', step-1, ' --> Step', step, 'Delta R-Sqrd p-val = ', delta_rsqr_p)
    
    print('-'*80) 