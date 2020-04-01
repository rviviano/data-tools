# Conditional Latent Growth Structural Equation Model Power Analysis

# Raymond Viviano
# rayviviano@gmail.com
# February 12th, 2020

# TODO: Vary models by small, medium, and large residual variance

suppressPackageStartupMessages(library(lavaan))
suppressPackageStartupMessages(library(mnormt))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(RColorBrewer))

options(warn=-1)

################### CHANGE MODEL AND SCRIPT PARAMETERS HERE ####################

intercept_mean_control <- 0 # Doesn't currently do anything.
intercept_mean_mci     <- 0

slope_mean_control <- 0
max_slope_mean_mci <- 1

intercept_variance <- 1
max_slope_variance <- 1

intercept_slope_covariance <- .2

residual_variances <- .5

intercept_age_coef_unstd    <-  .2
intercept_edu_coef_unstd    <- -.2
intercept_sex_coef_unstd    <-  .1
intercept_accult_coef_unstd <-  .1

slope_age_coef_unstd    <-  .2
slope_edu_coef_unstd    <- -.2
slope_sex_coef_unstd    <-  .1
slope_accult_coef_unstd <-  .1

nsimulations = 5000

##################### END USE DEFINED PARAMETERS SECTION  ######################


# Set Base Model, i.e., componenets that do not vary between models
growth_model_base <- "i =~ 1*y1 + 1*y2 + 1*y3\ns =~ 0*y1 + 1*y2 + 2*y3\n"
growth_model_base <- paste0(growth_model_base, "i ~~ ", intercept_variance, "*i\n")


# Age, edu, sex, and acculturation path coefficients
age_edu_terms <- paste0("i ~ ", intercept_age_coef_unstd, "*age\n")
age_edu_terms <- paste0(age_edu_terms, "s ~ ", slope_age_coef_unstd, "*age\n")
age_edu_terms <- paste0(age_edu_terms, "i ~ ", intercept_edu_coef_unstd, "*edu\n")
age_edu_terms <- paste0(age_edu_terms, "s ~ ", slope_edu_coef_unstd, "*edu\n")
age_edu_terms <- paste0(age_edu_terms, "i ~ ", intercept_sex_coef_unstd, "*sex\n")
age_edu_terms <- paste0(age_edu_terms, "s ~ ", slope_sex_coef_unstd, "*sex\n")
age_edu_terms <- paste0(age_edu_terms, "i ~ ", intercept_accult_coef_unstd, "*acculturation\n")
age_edu_terms <- paste0(age_edu_terms, "s ~ ", slope_accult_coef_unstd, "*acculturation\n")

# Set slope variance/covariance and residual variance for Ys
for(i in c(max_slope_variance/2, max_slope_variance)){
    # Vary slope variance between .1 and .5
    eval(parse(text=paste0('lgm_svar', i, ' <-paste0(growth_model_base, ',
                           '"s ~~ ', i, '*s +', intercept_slope_covariance, 
                           '*i\n")')))
    
    # Set residual variances to .5
    eval(parse(text=paste0('lgm_svar', i, ' <-paste0(lgm_svar', i, 
                           ', "\ny1 ~~ ', residual_variances,'*y1\ny2 ~~ ', 
                           residual_variances,'*y2\n",', '"y3 ~~ ', 
                           residual_variances,'*y3\n")')))  

    # Add age and edu covariates
    eval(parse(text=paste0('lgm_svar', i, ' <- paste0(lgm_svar', i, 
                           ', "\n", age_edu_terms)')))
}

# Set mean slope for MCI group
for(i in c(max_slope_variance/2, max_slope_variance)){
    for(j in seq(0, max_slope_mean_mci, by=max_slope_mean_mci/10)){
        # Add path coefficients between mci and slope (vary slope means
        # for the mci group across simulations)
        eval(parse(text=paste0('lgm_svar', i, '_mci_slp_mean', j,
                               ' <- paste0(lgm_svar', i, 
                               ', "\ni~0*1\ns~', j, '*1\n")')))
    }
}

# Test Model - Correctly Specified
m1 <- "
i =~ 1*y1 + 1*y2 + 1*y3
s =~ 0*y1 + 1*y2 + 2*y3
i ~~ i + s
s ~~ s
i ~ mci + age + edu + sex + acculturation
s ~ mci + age + edu + sex + acculturation
"

# Test Model - Misspecification 1 - Only MCI
m2 <- "
i =~ 1*y1 + 1*y2 + 1*y3
s =~ 0*y1 + 1*y2 + 2*y3
i ~~ i + s
s ~~ s
i ~ mci 
s ~ mci
"

# Test Model - Misspecification 2 - Add test site and language
m3 <- "
i =~ 1*y1 + 1*y2 + 1*y3
s =~ 0*y1 + 1*y2 + 2*y3
i ~~ i + s
s ~~ s
i ~ mci + age + edu + sex + site + acculturation + language
s ~ mci + age + edu + sex + site + acculturation + language
"

# Create dataframes to store simulation results
sim_df_m1 <- data.frame(matrix(ncol = 5, nrow = 0))
sim_df_m2 <- data.frame(matrix(ncol = 5, nrow = 0))
sim_df_m3 <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(sim_df_m1) <- c('mciSlopeMean', 's_mci', 'svar', 'r05', 'r05p05')
colnames(sim_df_m2) <- c('mciSlopeMean', 's_mci', 'svar', 'r05', 'r05p05')
colnames(sim_df_m3) <- c('mciSlopeMean', 's_mci', 'svar', 'r05', 'r05p05')

# Run the simulations
for(i in c(max_slope_variance/2, max_slope_variance)){
    for(j in seq(0, max_slope_mean_mci, by=max_slope_mean_mci/10)){
        cat(i,j,'\n')
        eval(parse(text=paste0('lgm_con <- lgm_svar', i, '_mci_slp_mean', j*0.0))) 
        eval(parse(text=paste0('lgm_mci <- lgm_svar', i, '_mci_slp_mean', j))) 
        
        rmsea05_count  <- 0
        r05_path05 <- 0
        mci_slope_path <- 0

        for(k in 1:nsimulations){
            # Simulate 100 MCI and 100 Control Observations
            df_con <- lavaan::simulateData(model=lgm_con, model.type="growth",
                                           sample.nobs=100)

            df_mci <- lavaan::simulateData(model=lgm_mci, model.type="growth",
                                           sample.nobs=100)

            # Add MCI columns to both data frames
            df_con[,"mci"] <- 0
            df_mci[,"mci"] <- 1

            # Combine Control and MCI datasets
            df <- rbind(df_con, df_mci)

            # Test the correctly specified and misspecified models
            m1_fit <- growth(m1, data=df, missing='fiml') 
            # m2_fit <- growth(m2, data=df, missing='fiml')
            # m3_fit <- growth(m3, data=df, missing='fiml')

            # Get RMSEA, s~mci estimate, and s~mci pvalue - Model 1
            smci_est <- as.double(parameterEstimates(m1_fit, standardized = FALSE)[15,'est'])
            pval <- as.double(parameterEstimates(m1_fit, standardized = FALSE)[15,'pvalue'])
            rmsea <- as.double(fitMeasures(m1_fit, c('rmsea'))[1])

            mci_slope_path <- mci_slope_path + smci_est

            if(rmsea <= .05){
                rmsea05_count <- rmsea05_count + 1
                if(pval <= .05){
                    r05_path05 <- r05_path05 + 1
                }
            } 
        }

        rmsea05_prop <- rmsea05_count/nsimulations
        r05_path05_prop <- r05_path05/nsimulations
        mci_slope_path_avg <- mci_slope_path/nsimulations

        temp_df <- data.frame(j, mci_slope_path_avg, i, rmsea05_prop, r05_path05_prop)
        colnames(temp_df) <- c('mciSlopeMean', 's_mci', 'svar', 'r05', 'r05p05')

        sim_df_m1 <- rbind(sim_df_m1, temp_df)
            
    }
}

# Sort data
sim_df_m1 <- sim_df_m1[order(sim_df_m1$mciSlopeMean, sim_df_m1$s_mci),]

print(sim_df_m1)

save(sim_df_m1, file="./sim_df_m1.Rda")