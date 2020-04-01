# Parallel Growth Power Analysis
# Raymond Viviano
# rayviviano@gmail.com
# February 16th, 2020

suppressPackageStartupMessages(library(lavaan))
suppressPackageStartupMessages(library(mnormt))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(RColorBrewer))

options(warn=-1)

################### CHANGE MODEL AND SCRIPT PARAMETERS HERE ####################
# Measurement model and covariance paths for subjective cognitive ability
scd_intercept_mean_control <- 0 
scd_intercept_mean_mci <- 0
scd_slope_mean_control <- 0
scd_slope_mean_mci <- 2 # Vary between 1 and 2
scd_intercept_variance <- 1
scd_slope_variance <- 1
scd_residual_variance <- .5
scd_intercept_slope_cov <- .2
scd_intercept_age_coef_unstd <-  .2
scd_intercept_edu_coef_unstd <- -.2
scd_intercept_sex_coef_unstd <-  .1
scd_intercept_accult_coef_unstd <-  .1
scd_slope_age_coef_unstd <-  .2
scd_slope_edu_coef_unstd <- -.2
scd_slope_sex_coef_unstd <-  .1
scd_slope_accult_coef_unstd <- .1

# Measurement model and covariance paths for objective cognitive performance
obj_intercept_mean_control <- 0 
obj_intercept_mean_mci <- 0
obj_slope_mean_control <- 0
# obj_slope_mean_mci varied below in the for loops
obj_intercept_variance <- 1
obj_slope_variance <- 1
obj_residual_variance <- .5
obj_intercept_slope_cov <- .2
obj_intercept_age_coef_unstd <- -.2
obj_intercept_edu_coef_unstd <- .1
obj_intercept_sex_coef_unstd <- .1
obj_intercept_accult_coef_unstd <- .1
obj_slope_age_coef_unstd <- -.2
obj_slope_edu_coef_unstd <- .1
obj_slope_sex_coef_unstd <- .1
obj_slope_accult_coef_unstd <- .1

# Structural model linking growth terms between subjective and objective cog
scd_int_to_obj_int <- -.3 
# scd_int_to_obj_slp varied in the for loops
# scd_slp_to_obj_slp varied in the for loops

# Other - General Simulation Parameters
nsimulations <- 5000

#################### END USER DEFINED PARAMETERS SECTION  ######################

# # Variables that the script varies...
# obj_slope_mean_mci; varied between -1 and -2
# scd_int_to_obj_slp; varied between .0 and -1
# scd_slp_to_obj_slp; variedbetween  .0 and -1

sim_df <- data.frame(matrix(ncol = 7, nrow = 0))

colnames(sim_df) <- c("scd_slope_mean_mci", "obj_slope_mean_mci", "scd_int_to_obj_int",
                      "scd_int_to_obj_slp", "scd_slp_to_obj_slp", "rmsea05_si_path05_prop",
                      "rmsea05_ss_path05_prop")

for (obj_slope_mean_mci in c(-1,-2)){
    for (scd_int_to_obj_slp in c(0, -.2, -.4, -.6, -.8, -1)){
        for (scd_slp_to_obj_slp in c(0, -.2, -.4, -.6, -.8, -1)){
            # Measurment model, latent means, vars, and covars for scd latent growth process
            # Reminder, latent slope varies by control/MCI status, so define that later
            growth_model_base <- "i_scd =~ 1*y1 + 1*y2 + 1*y3\n"
            growth_model_base <- paste0(growth_model_base, "s_scd =~ 0*y1 + 1*y2 + 2*y3\n")
            growth_model_base <- paste0(growth_model_base, "i_scd ~~ ", scd_intercept_variance, "*i_scd\n")
            growth_model_base <- paste0(growth_model_base, "s_scd ~~ ", scd_slope_variance, "*s_scd\n")
            growth_model_base <- paste0(growth_model_base, "i_scd ~~ ", scd_intercept_slope_cov, "*s_scd\n")

            # Measurment model, latent means, vars, and covars for objective cognitive
            # performance latent growth process. Reminder, latent slope varies by control 
            # or MCI status, so define that later
            growth_model_base <- paste0(growth_model_base, "i_obj =~ 1*x1 + 1*x2 + 1*x3\n")
            growth_model_base <- paste0(growth_model_base, "s_obj =~ 0*x1 + 1*x2 + 2*x3\n")
            growth_model_base <- paste0(growth_model_base, "i_obj ~~ ", obj_intercept_variance, "*i_obj\n")
            growth_model_base <- paste0(growth_model_base, "s_obj ~~ ", obj_slope_variance, "*s_obj\n")
            growth_model_base <- paste0(growth_model_base, "i_obj ~~ ", obj_intercept_slope_cov, "*s_obj\n")

            # Set residual variances for repeated measures of scd
            growth_model_base <- paste0(growth_model_base, "y1 ~~ ", scd_residual_variance, "*y1\n")
            growth_model_base <- paste0(growth_model_base, "y2 ~~ ", scd_residual_variance, "*y2\n")
            growth_model_base <- paste0(growth_model_base, "y3 ~~ ", scd_residual_variance, "*y3\n")

            # Set residual variances for repeated measures of objective performance
            growth_model_base <- paste0(growth_model_base, "x1 ~~ ", obj_residual_variance, "*x1\n")
            growth_model_base <- paste0(growth_model_base, "x2 ~~ ", obj_residual_variance, "*x2\n")
            growth_model_base <- paste0(growth_model_base, "x3 ~~ ", obj_residual_variance, "*x3\n")

            # Age, edu, sex, and acculturation path coefficients to scd latent processes
            growth_model_base <- paste0(growth_model_base, "i_scd ~ ", scd_intercept_age_coef_unstd, "*age\n")
            growth_model_base <- paste0(growth_model_base, "s_scd ~ ", scd_slope_age_coef_unstd, "*age\n")
            growth_model_base <- paste0(growth_model_base, "i_scd ~ ", scd_intercept_edu_coef_unstd, "*edu\n")
            growth_model_base <- paste0(growth_model_base, "s_scd ~ ", scd_slope_edu_coef_unstd, "*edu\n")
            growth_model_base <- paste0(growth_model_base, "i_scd ~ ", scd_intercept_sex_coef_unstd, "*sex\n")
            growth_model_base <- paste0(growth_model_base, "s_scd ~ ", scd_slope_sex_coef_unstd, "*sex\n")
            growth_model_base <- paste0(growth_model_base, "i_scd ~ ", scd_intercept_accult_coef_unstd, "*acculturation\n")
            growth_model_base <- paste0(growth_model_base, "s_scd ~ ", scd_slope_accult_coef_unstd, "*acculturation\n")

            # Age, edu, sex, and acculturation path coefficients to objective
            # cognitive performance latent processes
            growth_model_base <- paste0(growth_model_base, "i_obj ~ ", obj_intercept_age_coef_unstd, "*age\n")
            growth_model_base <- paste0(growth_model_base, "s_obj ~ ", obj_slope_age_coef_unstd, "*age\n")
            growth_model_base <- paste0(growth_model_base, "i_obj ~ ", obj_intercept_edu_coef_unstd, "*edu\n")
            growth_model_base <- paste0(growth_model_base, "s_obj ~ ", obj_slope_edu_coef_unstd, "*edu\n")
            growth_model_base <- paste0(growth_model_base, "i_obj ~ ", obj_intercept_sex_coef_unstd, "*sex\n")
            growth_model_base <- paste0(growth_model_base, "s_obj ~ ", obj_slope_sex_coef_unstd, "*sex\n")
            growth_model_base <- paste0(growth_model_base, "i_obj ~ ", obj_intercept_accult_coef_unstd, "*acculturation\n")
            growth_model_base <- paste0(growth_model_base, "s_obj ~ ", obj_slope_accult_coef_unstd, "*acculturation\n")

            # Create group specific models from the base model
            growth_model_con <- growth_model_base
            growth_model_mci <- growth_model_base

            # Set group specific paramters for different group models
            growth_model_con <- paste0(growth_model_con, "i_scd ~ ", scd_intercept_mean_control, "*1\n")
            growth_model_mci <- paste0(growth_model_mci, "i_scd ~ ", scd_intercept_mean_mci, "*1\n")
            growth_model_con <- paste0(growth_model_con, "s_scd ~ ", scd_slope_mean_control, "*1\n")
            growth_model_mci <- paste0(growth_model_mci, "s_scd ~ ", scd_slope_mean_mci, "*1\n")
            growth_model_con <- paste0(growth_model_con, "i_obj ~ ", obj_intercept_mean_control, "*1\n")
            growth_model_mci <- paste0(growth_model_mci, "i_obj ~ ", obj_intercept_mean_mci, "*1\n")
            growth_model_con <- paste0(growth_model_con, "s_obj ~ ", obj_slope_mean_control, "*1\n")
            growth_model_mci <- paste0(growth_model_mci, "s_obj ~ ", obj_slope_mean_mci, "*1\n")

            # Define structural model. Add paths between growth terms for con and mci models
            growth_model_con <- paste0(growth_model_con, "i_obj ~ ", scd_int_to_obj_int, "*i_scd\n")
            growth_model_mci <- paste0(growth_model_mci, "i_obj ~ ", scd_int_to_obj_int, "*i_scd\n")
            growth_model_con <- paste0(growth_model_con, "s_obj ~ ", scd_int_to_obj_slp, "*i_scd\n")
            growth_model_mci <- paste0(growth_model_mci, "s_obj ~ ", scd_int_to_obj_slp, "*i_scd\n")
            growth_model_con <- paste0(growth_model_con, "s_obj ~ ", scd_slp_to_obj_slp, "*s_scd\n")
            growth_model_mci <- paste0(growth_model_mci, "s_obj ~ ", scd_slp_to_obj_slp, "*s_scd\n")

            # Test Model - Correctly Specified
            m1 <- "
            i_scd =~ 1*y1 + 1*y2 + 1*y3
            s_scd =~ 0*y1 + 1*y2 + 2*y3
            i_scd ~~ i_scd + s_scd
            s_scd ~~ s_scd

            i_obj =~ 1*x1 + 1*x2 + 1*x3
            s_obj =~ 0*x1 + 1*x2 + 2*x3
            i_obj ~~ i_obj + s_obj
            s_obj ~~ s_obj

            i_scd ~ mci + age + edu + sex + acculturation
            s_scd ~ mci + age + edu + sex + acculturation
            i_obj ~ mci + age + edu + sex + acculturation + i_scd
            s_obj ~ mci + age + edu + sex + acculturation + i_scd + s_scd
            "

            # Run the simulations, SI <- Slope ON Intercept, SS <- Slope ON Slope
            rmsea05_si_path05 <- 0
            rmsea05_ss_path05 <- 0

            # Keep track of successful simulations without convergence issues
            ngoodsims <- 0

            for(k in 1:nsimulations){
                # Simulate 100 MCI and 100 Control Observations
                df_con <- lavaan::simulateData(model=growth_model_con, model.type="growth",
                                                sample.nobs=100)

                df_mci <- lavaan::simulateData(model=growth_model_mci, model.type="growth",
                                                sample.nobs=100)

                # Add MCI columns to both data frames
                df_con[,"mci"] <- 0
                df_mci[,"mci"] <- 1

                # Combine Control and MCI datasets
                df <- rbind(df_con, df_mci)

                # Test the correctly specified model
                e <- try(fit <- lavaan::growth(m1, data=df, missing='fiml'))

                # Skip this simulation if model did not converge
                if(inherits(e, "try-error")) {
                    fit <- NULL
                }

                # Check fit errors/model convergence errors
                e2 <- try(fitStatsDummy <- fitmeasures(fit))
                if(inherits(e2, "try-error")) {
                    fit <- NULL
                }

                # If the model ran successfully, extract fit stats and parameter estimates
                if (!is.null(fit)) {
                    # Increment number of successful simulations
                    ngoodsims <- ngoodsims + 1

                    # Get RMSEA, s_obj~i_scd and s_obj~s_scd pvalues
                    rmsea <- as.double(fitMeasures(fit, c('rmsea'))[1])
                    si_pval <- as.double(parameterEstimates(fit, standardized = FALSE)[40, 'pvalue'])
                    ss_pval <- as.double(parameterEstimates(fit, standardized = FALSE)[41, 'pvalue'])

                    if(rmsea <= .05){
                        if(si_pval <= .05){
                            rmsea05_si_path05 <- rmsea05_si_path05 + 1
                        } else if(ss_pval <= .05){
                        rmsea05_ss_path05 <- rmsea05_ss_path05 + 1 
                        }
                    }
                }
            }

            rmsea05_si_path05_prop <- rmsea05_si_path05/ngoodsims
            rmsea05_ss_path05_prop <- rmsea05_ss_path05/ngoodsims

            # Create dataframe to store simulation results
            temp <- data.frame(scd_slope_mean_mci, obj_slope_mean_mci, scd_int_to_obj_int,
                                    scd_int_to_obj_slp, scd_slp_to_obj_slp, rmsea05_si_path05_prop,
                                    rmsea05_ss_path05_prop)

            colnames(temp) <- c("scd_slope_mean_mci", "obj_slope_mean_mci", "scd_int_to_obj_int",
                                "scd_int_to_obj_slp", "scd_slp_to_obj_slp", "rmsea05_si_path05_prop",
                                "rmsea05_ss_path05_prop")

            sim_df <- rbind(sim_df, temp)

        }
    }
}

foo <- as.character(Sys.time())

print(foo)
print(paste0("./sim_parallel_growth-", as.character(Sys.time()),".Rda"))

save(sim_df, file=paste0("./sim_parallel_growth-", as.character(Sys.time()),".Rda"))
write.csv(sim_df, paste0("./sim_parallel_growth-", as.character(Sys.time()),".csv"))