# Binomial logistic regression with mediation simulation code
# NOTE: The assumption is that the mediator is a continuous variable

# Raymond Viviano
# rayviviano@gmail.com
# March 25, 2020

# TODO: Implement checks to make sure that all models converged


#' Inverse Logit Function
inv.logit <- function(p){
    return(exp(p)/(1+exp(p)))
}


#' Sobel indirect effect standard error approximation
sobel.std.err <- function(a, b, se.a, se.b){
    return((b^2 * se.a^2) + (a^2 + se.b^2))
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

# Some equations to keep in mind for logistic regression with mediation:
#
# EQ1.  Y   =   B0_1 + Tau*X + e1
# EQ2.  Y   =   B0_2 + Tau_Prime*X + Beta*M + e2
# EQ3.  M   =   B0_3 + alpha*X + e3
# EQ4.  Tau =   alpha*beta + Tau_Prime   || (c = ab + c') (approx for logistic)
#
# Tau = total effect, tau_prime = direct effect of X on Y accounting for
# mediating relationship, alpha = effect of X on M, beta = effect of M on Y
# B0_1,2,&3 are intercept terms for their respective equations. 
#
# Assuming zero error, the mean, or expected value, for the intercept term for 
# EQ1 can be defined by the other variables.
#
# EQ5. E[B0_1] = E[B0_2] + Tau_Prime*E[X] + Beta*(E[B0_3] + Alpha*E[X]) - Tau*[X]
#
# If we force the mean of X to 1, then we get equation 6:
#
# EQ6. E[B0_1] = E[B0_2] + Tau_Prime + Beta*E[B0_3] + Alpha*Beta - Tau
#
# As Tau_Prime - Tau = -1 * Alpha*Beta, we get equation 7:
#
# EQ7. E[B0_1] = E[B0_2] + Beta*E[B0_3]
#
# We also get equation 7 if we force the mean of x to 0...and the derivation is easier...
#
# Define B0_1, B0_2, B0_3, alpha, beta, tau, and tau_prime for the data
# generation process.
#
# Use Eq3 to generate M, then use Eq2 to generate Y
# Test that you recover the expected values for alpha, beta, tau, tau_prime,
# B0_1, B02, and B03

# Specify indirect and direct effect terms
alpha = .2
beta = .2
tau_prime = .1

# Specify total effect (note that for logit this is actually approximate)
tau = alpha*beta + tau_prime

# Intercept terms for EQ2 and EQ3 (Means or expected values)       
B0_2 <- .2        
B0_3 <- .4 

# Calculate intercept term for EQ1 based on the other specified values
B0_1 <- B0_2 + beta*B0_3

# Ensure reproducible results
set.seed(10000)

# Define sample size
n <- 16000

# Define number of simulations
nsims <- 50

# Define matrix to hold simulation data (estimates, std.errors, and pvals)
sim.prms <- matrix(NA, nrow=nsims, ncol=20)

# Loop through simulations
for(i in 1:nsims){
    # Generate independant variable data, uniform random between -1 and 1
    # The expected value, or expected mean of this vector should be 0
    x <- runif(n, -1, 1)

    # Generate data for the mediating variable, the mean or expected value 
    # should equal B03. But, we'll add some noise
    m <- alpha*x + B0_3 
    # Add noise
    m <- m + rnorm(length(m), 0, .2)

    # True data-generation-process for Bernoulli trials
    y <- rbinom(n, 1, inv.logit(B0_2 + tau_prime*x + beta*m))

    # Estimate logit model
    model1 <- glm(y ~ x, family=binomial(link=logit))
    model2 <- glm(y ~ x + m, family=binomial(link=logit))
    model3 <-  lm(m ~ x)

    # TODO: Implement checks to make sure that all models converged

    # Put model paramters in simulation matrix
    sim.prms[i, 1] <- model1$coef[1]                    # Estimate for B0_1
    sim.prms[i, 2] <- model2$coef[1]                    # Estimate for B0_2
    sim.prms[i, 3] <- model3$coef[1]                    # Estimate for B0_3
    sim.prms[i, 4] <- model3$coef[2]                    # Estimate for Alpha
    sim.prms[i, 5] <- model2$coef[3]                    # Estimate for Beta
    sim.prms[i, 6] <- sim.prms[i, 4] * sim.prms[i, 5]   # Estimate for Indirect 
    sim.prms[i, 7] <- model1$coef[2]                    # Estimate for Tau
    sim.prms[i, 8] <- model2$coef[2]                    # Estimate for Tau_Prime
    sim.prms[i, 9] <- summary(model3)$coefficients[,2][2]  # Alpha Std.Err
    sim.prms[i,10] <- summary(model2)$coefficients[,2][3]  # Beta Std.Err
    sim.prms[i,11] <- sobel.std.err(sim.prms[i,4], sim.prms[i,5], 
                                    sim.prms[i,9], sim.prms[i,10]) # Indirect Std.Err
    sim.prms[i,12] <- summary(model1)$coefficients[,2][2]  # Tau Std.Err
    sim.prms[i,13] <- summary(model2)$coefficients[,2][2]  # Tau_Prime Std.Err
    sim.prms[i,14] <- summary(model1)$coefficients[,4][1]  # P-Value for B0_1
    sim.prms[i,15] <- summary(model2)$coefficients[,4][1]  # P-Value for B0_2
    sim.prms[i,16] <- summary(model3)$coefficients[,4][1]  # P-Value for B0_3
    sim.prms[i,17] <- summary(model3)$coefficients[,4][2]  # P-Value for Alpha
    sim.prms[i,18] <- summary(model2)$coefficients[,4][3]  # P-Value for Beta
    sim.prms[i,19] <- summary(model1)$coefficients[,4][2]  # P-Value for Tau
    sim.prms[i,20] <- summary(model2)$coefficients[,4][2]  # P-Value for Tau'

}

# Note: The assumption for this simulation is that the mediator is continuous
cat("Note: This simulation assumes that the mediator is continous...\n\n")

# Get Covarage probabilities for the parameters
# B0_1.coverage <- coverage.prob(sim.prms[,1], sim.prms[,], B0_1, .95, n-model1$rank) TODO: Add these std.errs to sim.prms
# B0_2.coverage <- coverage.prob(sim.prms[,2], sim.prms[,], B0_2, .95, n-model2$rank) TODO: Add these std.errs to sim.prms
# B0_3.coverage <- coverage.prob(sim.prms[,3], sim.prms[,], B0_3, .95, n-model3$rank) TODO: Add these std.errs to sim.prms
alpha.coverage     <- coverage.prob(sim.prms[,4], sim.prms[, 9], alpha, .95, n-model3$rank)
beta.coverage      <- coverage.prob(sim.prms[,5], sim.prms[,10], beta,  .95, n-model2$rank)
indirect.coverage  <- coverage.prob(sim.prms[,6], sim.prms[,11], alpha*beta, .95, n-model2$rank) # TODO: Figure out confidence interval for sobel std.err
tau.coverage       <- coverage.prob(sim.prms[,7], sim.prms[,12], tau,        .95, n-model3$rank)
tau.prime.coverage <- coverage.prob(sim.prms[,8], sim.prms[,13], tau_prime,  .95, n-model2$rank)

# # Intercept coverage probabilities     -- TODO
# print("Coverage Probability for B0_1")
# print(B0_1.coverage)
# print("Coverage Probability for B0_2")
# print(B0_2.coverage)
# print("Coverage Probability for B0_3")
# print(B0_3.coverage)

print("alpha.coverage")
print(alpha.coverage)
print("beta.coverage")     
print(beta.coverage) 
print("indirect.coverage")      
print(indirect.coverage)  
print("tau.coverage")
print(tau.coverage)        
print("tau.prime.coverage")
print(tau.prime.coverage)  


# TODO: Update/Fix Proportion Code
# power.detect.effect(sim.prms[,5])
# power.detect.effect(sim.prms[,6])