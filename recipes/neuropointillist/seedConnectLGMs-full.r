suppressWarnings(suppressPackageStartupMessages(library(lavaan)))
suppressWarnings(suppressPackageStartupMessages(library(tidyr)))

processVoxel <- function(v) {
    # Get all participant connectivity values for voxel v
    Y <- voxeldat[,v]

    # Merge covariates with voxel data, designmat is global in the npoint script
    df <- data.frame(designmat$idnum, designmat$time, designmat$mfq, 
                     designmat$age, Y)

    # Rename cols to sensible variables. E.g., "designmat.idnum" becomes "idnum"
    colnames(df) <- c("idnum", "time", "mfq", "age", "Y")
    
    # Widen the data for lavaan, no time-varying covariates to worry about here
    df_w <- df[,] %>% 
        gather(variable, value, Y) %>%
        unite(var_time, variable, time) %>%
        spread(var_time, value)

    require(lavaan)

    # Specify the LGM
    cor_lgc_model <- '
    # Intercept and Slope, Time unit Approx 18 months
    i =~ 1*Y_1 + 1*Y_2 + 1*Y_3
    s =~ 0*Y_1 + 1*Y_2 + 2*Y_3

    # # Regressions
    i ~ age + mfq
    s ~ age + mfq

    # Variances and covariances explicitly coded
    i ~~ s 
    '

    # Try to fit the LGM for the voxel
    e <- try(fit <- lavaan::growth(cor_lgc_model, data=df_w, 
                                   estimator="MLR", missing="ML"))

    # One bad voxel should spoil the whole bunch. Return -999 for everything
    if(inherits(e, "try-error")) {
        message("error thrown at voxel",v)
        message(e)
        fit <- NULL
    }

    # Check residual variances for negative values. If the negative variances
    # are small enough though and SEs suggest that values could be greater than
    # 0, the negative variance could be due to sampling error - constrain the 
    # residual variance to 0 and try again
    if (!is.null(fit)) {
        unstd_ests <- parameterEstimates(fit)
        vars <- subset(unstd_ests, op=="~~" & lhs==rhs)
        vars <- na.omit(vars)
        # If there are small negative variances constrain to positive value 
        # and run again, else set fit to null and return -999 for everything
        if (sum(vars$est < 0) > 0) { 
            # Create an updated model
            updated_model <- cor_lgc_model
            # Run through all variances and covariances
            for (i in 1:nrow(vars)){
                # Check if latent intercept variance < 0
                if ((vars[i,]$est < 0)) {
                    # If variance negative and the 95% confidence interval 
                    # suggests that variance could be greater than 0, then
                    # constrain the residual variance to 0
                    if ((vars[i,]$est + vars[i,]$se*1.96 > 0)) {
                        # Add residual variance constraint to model
                        updated_model <- paste0(updated_model, "\n", 
                                                vars[i,]$lhs, "~~0*", 
                                                vars[i,]$rhs)
                    }
                } 
            }     
  
            e <- try(fit <- lavaan::growth(updated_model, data=df_w, 
                                            estimator="MLR", missing="ML"))

            if(inherits(e, "try-error")) {
                message("error thrown at voxel",v)
                message(e)
                fit <- NULL
            } 
        }
    }

    # If there are still negative latent variances. The voxel is beyond saving
    if (!is.null(fit)) {
        unstd_ests <- parameterEstimates(fit)
        vars <- subset(unstd_ests, op=="~~" & lhs==rhs)
        if (sum(vars$est < 0) > 0) { 
            # Negative variances Exist. Set fit to NULL and return -999 for all
            fit <- NULL
        }
    }

    # Return -999 if model did not converge
    e2 <- try(fitStatsDummy <- fitmeasures(fit))
    if(inherits(e2, "try-error")) {
        message("convergence error thrown at voxel", v)
        message(e2)
        fit <- NULL
    }

    # If the model ran successfully, extract fit stats and parameter estimates
    if (!is.null(fit)) {
        # Unstandardized Parameter Estimates
        unstd_ests <- parameterEstimates(fit)

        # Standardized Parameter Estimates
        param_ests <- standardizedsolution(fit)

        # Name parameters/pathways in the measurement and structural model
        unstd_ests$name <- paste0(unstd_ests$lhs, unstd_ests$op, unstd_ests$rhs)
        param_ests$name <- paste0(param_ests$lhs, param_ests$op, param_ests$rhs)

        # Get fit measures of interest
        fitStats <- fitmeasures(fit, c("chisq", "df.scaled", "pvalue", "rmsea", 
                                       "rmsea.ci.lower", "rmsea.ci.upper", 
                                       "cfi", "tli", "srmr"))

        # Return unstandardized and standardized estimates 
        latent_int <- unstd_ests[unstd_ests$lhs %in% c('i', 's') & 
                                 unstd_ests$op %in% c('~1'),]
        latent_reg <- unstd_ests[unstd_ests$lhs %in% c('i', 's') & 
                                 unstd_ests$op %in% c('~'),]
        latent_varcov <- unstd_ests[unstd_ests$lhs %in% c('i', 's') & 
                                    unstd_ests$op %in% c('~~'),]

        latent_int_std <- param_ests[param_ests$lhs %in% c('i', 's') & 
                                     param_ests$op %in% c('~1'),]
        latent_reg_std <- param_ests[param_ests$lhs %in% c('i', 's') & 
                                     param_ests$op %in% c('~'),]
        latent_varcov_std <- param_ests[param_ests$lhs %in% c('i', 's') & 
                                        param_ests$op %in% c('~~'),]

        # Bind all coefficients of interest together
        ret <- rbind(latent_int, latent_reg, latent_varcov)
        ret_std <- rbind(latent_int_std, latent_reg_std, latent_varcov_std)

        # Make sure order of coefficients is correct
        ret <- ret[order(ret$name),]
        ret_std <- ret_std[order(ret_std$name),]

        # Create a new container for values to return to npoint, only keep
        # estimates and pvalues
        retvals <- c(ret$est, ret$pvalue, ret_std$est.std, ret_std$pvalue, fitStats)

        # Rename  the return values
        names(retvals) <-c(paste(ret$name, ".est", sep=""), 
                           paste(ret$name, ".pvalue", sep=""),
                           paste(ret_std$name,".est_std", sep=""), 
                           paste(ret_std$name, ".pvalue_std", sep=""), 
                           "chisq", "df", "pvalue", "rmsea", "rmsea.ci.lower", 
                           "rmsea.ci.upper", "cfi", "tli", "srmr")

    } else {
        # If the LGM didn't fit, e.g., from negative latent variances or from
        # other potential convergence issues, then return -999 for everything
        retvals <- rep(-999, 45)
        names(retvals) <- c("i~~i.est", "i~~s.est", "i~1.est", "i~age.est", 
                            "i~mfq.est", "s~~s.est", "s~1.est", "s~age.est", 
                            "s~mfq.est", "i~~i.pvalue", "i~~s.pvalue", 
                            "i~1.pvalue", "i~age.pvalue", "i~mfq.pvalue", 
                            "s~~s.pvalue", "s~1.pvalue", "s~age.pvalue", 
                            "s~mfq.pvalue", "i~~i.est_std", "i~~s.est_std", 
                            "i~1.est_std", "i~age.est_std", "i~mfq.est_std", 
                            "s~~s.est_std", "s~1.est_std", "s~age.est_std", 
                            "s~mfq.est_std", "i~~i.pvalue_std", 
                            "i~~s.pvalue_std", "i~1.pvalue_std", 
                            "i~age.pvalue_std", "i~mfq.pvalue_std",
                            "s~~s.pvalue_std", "s~1.pvalue_std", 
                            "s~age.pvalue_std", "s~mfq.pvalue_std", 
                            "chisq", "df", "pvalue", "rmsea", "rmsea.ci.lower", 
                            "rmsea.ci.upper", "cfi", "tli", "srmr")     
    }

    # Provide the return values to the main npoint script
    retvals
}
