# Multinomial logistic regression simulation code
# For K levels of the dependent variable, estimates K-1 sets of 
# coefficients for the independent variables.

# Currently set up for 3 classes but can be extended to any number

# Raymond Viviano
# rayviviano@gmail.com
# March 23, 2020

# TODO: Encapsulate certain code sections into smaller functions
# TODO: Bound coverage probability between 0 & 1

# Easy access to mlogit models
suppressPackageStartupMessages(library(mlogit))

# Set error reporting options
options(show.error.locations=TRUE)
options(error=recover)


#' Coverage Probability, i.e., the probability that the confidence interval (ci)
#' for a parameter will contain the true value
#' Inputs:
#'          b       - Estimate
#'          se      - Standard Error of the Estimate
#'          true    - True value of parameter for data-generating process
#'          cl      - Confidence level
#'          dof     - Degrees of Freedom
#' Output:
#'          List containing coverage probability, vect
coverage.prob <- function(b, se, true, cl=.95, dof=Inf){
    # Compute quantile based on confidence level
    coverage.quantile <- cl + (1-cl)/2

    # Compute confidence interval upper and lower bounds
    ci.lower <- b - qt(coverage.quantile, df=dof)*se
    ci.upper <- b + qt(coverage.quantile, df=dof)*se

    # For each ci generated for each b/se pair, eval if true param is in ci
    ci.contain.true <- ifelse(true>=ci.lower & true<=ci.upper, 1, 0)

    # Calculate coverage probability
    cp <- mean(ci.contain.true)

    # Calculate Monte Carlo Error
    mc.err.lower <- cp - 1.96*sqrt((cp*(1-cp))/length(b))
    mc.err.upper <- cp + 1.96*sqrt((cp*(1-cp))/length(b))

    # Return coverage probability and error
    return(list(cp=cp, ci=cbind(mc.err.lower, mc.err.upper)))
}


#' Power to detect any effect at various alpha levels
#' Takes a vector of p-values and calculates the proportion of modell where the 
#' p-val for the effect of interest was < .5, .01, and .001.
power.detect.effect <- function(p.val.vec){
    less.than.05 <- ifelse(p.val.vec < .05, 1, 0)
    prop.05 <- sum(less.than.05)/length(less.than.05)

    less.than.01 <- ifelse(p.val.vec < .01, 1, 0)
    prop.01 <- sum(less.than.01)/length(less.than.01)

    less.than.001 <- ifelse(p.val.vec < .001, 1, 0)
    prop.001 <- sum(less.than.001)/length(less.than.001)

    cat(paste0('Power to detect effect at .05: ', prop.05, '\n'))
    cat(paste0('Power to detect effect at .01: ', prop.01, '\n'))
    cat(paste0('Power to detect effect at .001: ', prop.001, '\n'))
}


# Ensure reproducible results
set.seed(10000)

# Set logistic regression parameters (Can be calculated as log(oddsRatio)) for
# Classes A and B (against reference C).

# Class a
b0a <- .2        # Set Logistic Regression Intercept True Value
b1a <- .5        # Set Logistic Regression Slope True Value 

# Class b
b0b <- -.2
b1b <- .75

# Define sample size
n <- 16000

# Define number of simulations
nsims <- 55

# Define matrix to hold simulation data (estimates, std.errors, and pvals)
sim.prms <- matrix(NA, nrow=nsims, ncol=12)

# Loop through simulations
for(i in 1:nsims){
    # Participant/subject/case/index
    case = c(1:n)

    # Generate independant variable data, uniform random between -1 and 1
    x <- runif(n, -1, 1)

    # Compute probabilities of each class based on generated data
    prob.a <- exp(b0a +b1a*x) / (1 + exp(b0a + b1a*x) + exp(b0b + b1b*x))
    prob.b <- exp(b0b +b1b*x) / (1 + exp(b0a + b1a*x) + exp(b0b + b1b*x))
    prob.c <- 1 - prob.a - prob.b

    # Generate y data based on calculated probabilities
    y <- rep(NA, n)
    for (j in 1:n){
        y[j] <- sample(c("a", "b", "c"), 1, replace = TRUE, 
                       prob = c(prob.a[j], prob.b[j], prob.c[j]))
    }

    # Prep data for mlogit package
    sim.data <- data.frame(case, x, y, stringsAsFactors=FALSE)

    # Convert data to long format for mlogit
    sim.data <- mlogit.data(sim.data, shape='wide', choice="y", id.var='case')

    # Estimate logit model
    model <- mlogit(y ~ -1|x, data=sim.data, reflevel='c')  

    # Get variance-covariance matrix
    var.covar <- vcov(model)

    # Put model paramters in simulation matrix
    sim.prms[i,1] <- model$coef[1]             # Estimate for b0a
    sim.prms[i,2] <- model$coef[3]             # Estimate for b1a
    sim.prms[i,3] <- model$coef[2]             # Estimate for b0b
    sim.prms[i,4] <- model$coef[4]             # Estimate for b1b
    sim.prms[i,5] <- sqrt(diag(var.covar)[1])  # b0a standard error
    sim.prms[i,6] <- sqrt(diag(var.covar)[3])  # b1a standard error
    sim.prms[i,7] <- sqrt(diag(var.covar)[2])  # b0b standard error
    sim.prms[i,8] <- sqrt(diag(var.covar)[4])  # b1b standard error
    sim.prms[i,9]  <- summary(model)$CoefTable[,4][1] # b0a pvalue
    sim.prms[i,10] <- summary(model)$CoefTable[,4][3] # b1a pvalue
    sim.prms[i,11] <- summary(model)$CoefTable[,4][2] # b0b pvalue
    sim.prms[i,12] <- summary(model)$CoefTable[,4][4] # b1b pvalue
}

# Get Covarage probabilities for the parameters
b0a.coverage <- coverage.prob(sim.prms[,1], sim.prms[,5], b0a, .95, n-length(coef(model)))
b1a.coverage <- coverage.prob(sim.prms[,2], sim.prms[,6], b1a, .95, n-length(coef(model)))
b0b.coverage <- coverage.prob(sim.prms[,3], sim.prms[,7], b0b, .95, n-length(coef(model)))
b1b.coverage <- coverage.prob(sim.prms[,4], sim.prms[,8], b1b, .95, n-length(coef(model)))

# Print intercept for class a coverage probability and 95% confidence interval
print("Coverage Probability for b0a")
print(b0a.coverage)

# Print coefficent on x for class a coverage probability and 95% confidence interval
print("Coverage Probability for b1a")
print(b1a.coverage)

# Print intercept for class b coverage probability and 95% confidence interval
print("Coverage Probability for b0b")
print(b0b.coverage)

# Print coefficent on x for class b coverage probability and 95% confidence interval
print("Coverage Probability for b1b")
print(b1b.coverage)

cat("Power to detect b0a\n")
power.detect.effect(sim.prms[,9])

cat("\nPower to detect b1a\n")
power.detect.effect(sim.prms[,10])

cat("Power to detect b0b\n")
power.detect.effect(sim.prms[,11])

cat("\nPower to detect b1b\n")
power.detect.effect(sim.prms[,12])