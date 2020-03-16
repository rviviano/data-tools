
# Get hazard ratios comparing different classes\groups. Code for both simple and 
# complex-survey data

# Raymond Viviano
# August 28th, 2019

library(foreign)
library(survey)
library(survival)


coxRegression <- function(df){
    #' Inputs:  dataframe of interest
    #'          integer representing class to analyze

    sDesign <- svydesign(id=~df$secu, strata=df$stratum, weights=df$gwgtr, 
                         data=df, nest=TRUE)
        
    sDesignSub <- subset(sDesign, df$subpop65==1)

    svycoxph(Surv(time, status) ~ Class, design=sDesign, 
             subset=subpop==1, data=df)
}

coxRegressionSurvey <- function(df){
    #' Inputs:  dataframe of interest
    #'          integer representing class to analyze

    sDesign <- svydesign(id=~df$secu, strata=df$stratum, weights=df$gwgtr, 
                         data=df, nest=TRUE)
        
    sDesignSub <- subset(sDesign, df$subpop65==1)

    svycoxph(Surv(time, status) ~ Class, design=sDesign, 
             subset=subpop==1, data=df)
}