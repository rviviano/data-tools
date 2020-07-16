# Binomial logistic regression simulation code

# Raymond Viviano
# rayviviano@gmail.com
# March 23, 2020


#' Inverse Logit Function
inv.logit <- function(p){
    return(exp(p)/(1+exp(p)))
}


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

# Set logistic regression parameters (Can be calculated as log(oddsRatio))
b0 <- .2        # Set Logistic Regression Intercept True Value
b1 <- .5        # Set Logistic Regression Slope True Value 

# Define sample size
n <- 16000

# Define number of simulations
nsims <- 20

# Define matrix to hold simulation data (estimates, std.errors, and pvals)
sim.prms <- matrix(NA, nrow=nsims, ncol=6)

# Loop through simulations
for(i in 1:nsims){
    # Generate independant variable data, uniform random between -1 and 1
    x <- runif(n, -1, 1)

    # Generate true data-generation-process Bernoulli trials
    y <- rbinom(n, 1, inv.logit(b0 + b1*x))

    # Estimate logit model
    model <- glm(y~x, family=binomial(link=logit))

    # Get variance-covariance matrix
    var.covar <- vcov(model)

    # Put model paramters in simulation matrix
    sim.prms[i,1] <- model$coef[1]             # Estimate for b0
    sim.prms[i,2] <- model$coef[2]             # Estimate for b1
    sim.prms[i,3] <- sqrt(diag(var.covar)[1])  # b0 standard error
    sim.prms[i,4] <- sqrt(diag(var.covar)[2])  # b1 standard error
    sim.prms[i,5] <- summary(model)$coefficients[,4][1] # b0 pvalue
    sim.prms[i,6] <- summary(model)$coefficients[,4][2] # b1 pvalue
}

# Get Covarage probabilities for the parameters
b0.coverage <- coverage.prob(sim.prms[,1], sim.prms[,3], b0, .95, n-model$rank)
b1.coverage <- coverage.prob(sim.prms[,2], sim.prms[,4], b1, .95, n-model$rank)

# Print intercept coverage probability and 95% confidence interval
print("Coverage Probability for b0")
print(b0.coverage)

# Print coefficent on x coverage probability and 95% confidence interval
print("Coverage Probability for b1")
print(b1.coverage)

# Power for intercept
cat("Power to detect b0\n")
power.detect.effect(sim.prms[,5])

cat("\nPower to detect b1\n")
# Power for coefficient on x
power.detect.effect(sim.prms[,6])