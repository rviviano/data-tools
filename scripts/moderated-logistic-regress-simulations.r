# Binomial logistic regression with moderation simulation code

# Raymond Viviano
# rayviviano@gmail.com
# March 27th, 2020

# TODO: Implement checks to make sure that all models converge

library(ggplot2)

########################### User Definition Section ############################

# Specify Intercept and Beta Terms 
B0 <- logit(0.1) # Intercept. 
B1 <- .2         # Variable 1 Beta
B2 <- .1         # Variable 2 Beta
B3.Odds.Seq <- seq(1.05, 2.55, by=.1) # Loop through ORs for interaction term

# Ensure reproducible results
set.seed(10000)

# Define sample size
n <- 16000

# Define number of simulations per B3 Odds Ratio
nsims <- 100

######################## End of User Definition Section ########################

############################# Function Definitions #############################

#' Logit function
logit <- function(p){
    return(log(as.double(p)/(1.0-as.double(p))))
}


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
#' Takes a vector of p-values and calculates the proportion of models where the 
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


#' Power to detect effect at .05 level based on simulations
power.at.05 <- function(p.val.vec){
    less.than.05 <- ifelse(p.val.vec < .05, 1, 0)
    return(sum(less.than.05)/length(less.than.05))
}


#' plotLines
#' 
#' Lineplots comparing continuous variable scores across levels of another var
#' 
#' Inputs: 
#'     df:       Dataframe with vars of interest.
#'     title:    Title above chart 
#'     outDir:   Directory to save plot to
plotLines <- function(df, title, outDir){
    # Create plot
    p <- ggplot(data=df, aes(x=B3.Odds.Seq, y=B3.Powr.Seq)) 

    # Add lines for each group
    p <- p + geom_line(stat='identity', size=.85, alpha=0.9) 
    
    # Unique shape for each group placed at each data point
    p <- p + geom_point(size=2.5, alpha=.45) 

    # Pad the distance between the x-axis and smallest plotted value slightly
    p <- p + scale_y_continuous(limits=c(0,1))

    # Pad/add extra space to the ends of the x-axis
    p <- p + scale_x_continuous(expand=c(0.05, 0.05)) 

    # Remove gray background, format x-axis text, format legend 
    p <- p + theme(panel.grid.major=element_blank(), 
                   panel.grid.minor=element_blank(), 
                   panel.background=element_blank(), 
                   axis.line=element_line(colour="black"), 
                   axis.text.x=element_text(angle=0, size=14),
                   axis.text.y=element_text(angle=0, size=14), 
                   legend.key=element_rect(color=NA, fill=NA),
                   legend.background=element_blank(),
                   legend.key.width=unit(3, "line"),
                   legend.justification=c(1, 1),
                   legend.position=c(.2, 1)) 
    
    # Add axis labels
    p <- p + labs(y='Proportion of Simulations with Significant Interaction Term at p < .05', 
                  x='Odds Ratio of the Interaction')  

    # Remove any components to the legend that could relate to alpha or size
    p <- p + scale_alpha(guide='none') 
    p <- p + scale_size(guide='none') 

    # Add a title to the plot
    p <- p + ggtitle(title) 

    # Save Plots
    pngOut <- paste(outDir, "/moderatedLogitInteractionPower.png", sep="")    
    svgOut <- paste(outDir, "/moderatedLogitInteractionPower.svg", sep="")   
    
    ggsave(pngOut, p)
    ggsave(svgOut, p)
}

# Check power at varying values for B3 (Interaction Beta)
B3.Powr.Seq <- rep(0, length(B3.Odds.Seq)) # Hold power for ORs, init all to 0

for(j in 1:length(B3.Odds.Seq)){
    # Take the natural log of the expected odds ratio to get interaction beta
    B3.Odds.Ratio <- B3.Odds.Seq[j]
    B3 <- log(B3.Odds.Ratio)

    # Define matrix to hold simulation data (estimates, std.errors, and pvals)
    sim.prms <- matrix(NA, nrow=nsims, ncol=12)

    # Print statements to track script progress
    cat(paste0("Calculating power for B3 = ", round(B3,2), ", exp(", round(B3,2), 
               ") = Odds Ratio of ", B3.Odds.Ratio, "\n\nIteration: "))

    # Loop through simulations
    for(i in 1:nsims){
        # Print statements to observe progress
        if ((i-1)%%10==0){
            cat(paste0(i,", "))
        }

        # Generate independant variable data, uniform random between -1 and 1
        x <- runif(n, -1, 1)

        # Generate normal dist data for the moderating variable
        m <- rnorm(n, 0, 1)

        # Element-wise multiply x and m for interaction
        xm <- x*m

        # Data-generation-process for Bernoulli trials
        y <- rbinom(n, 1, inv.logit(B0 + B1*x + B2*m + B3*xm))

        # Estimate logit model
        model <- glm(y ~ x + m + xm, family=binomial(link=logit))

        # TODO: Implement checks to make sure that models converge
        
        # Put model paramters in simulation matrix
        sim.prms[i, 1] <- model$coef[1]                       # Estimate for B0
        sim.prms[i, 2] <- model$coef[2]                       # Estimate for B1
        sim.prms[i, 3] <- model$coef[3]                       # Estimate for B2
        sim.prms[i, 4] <- model$coef[4]                       # Estimate for B3
        sim.prms[i, 5] <- summary(model)$coefficients[,2][1]  # StdError for B0
        sim.prms[i, 6] <- summary(model)$coefficients[,2][2]  # StdError for B1
        sim.prms[i, 7] <- summary(model)$coefficients[,2][3]  # StdError for B2
        sim.prms[i, 8] <- summary(model)$coefficients[,2][4]  # StdError for B3
        sim.prms[i, 9] <- summary(model)$coefficients[,4][1]  # P-Value for B0
        sim.prms[i,10] <- summary(model)$coefficients[,4][2]  # P-Value for B1
        sim.prms[i,11] <- summary(model)$coefficients[,4][3]  # P-Value for B2
        sim.prms[i,12] <- summary(model)$coefficients[,4][4]  # P-Value for B3

    }

    # Newline for prettier prints
    cat("\n")

    # Get Covarage probabilities for the interaction
    B3.coverage <- coverage.prob(sim.prms[,4], sim.prms[,8], B3, .95, n-model$rank)

    # Intercept coverage probabilities
    cat("Coverage Probability for B3 (Interaction Beta)\n")
    cat(paste0(B3.coverage$cp, " [", B3.coverage$ci[1,1], ",", 
               B3.coverage$ci[1,2], "]\n"))

    # Power to detect moderation
    print("Power to detect interaction term:")
    power.detect.effect(sim.prms[,12])

    # Update B3.Powr.Seq
    B3.Powr.Seq[j] <- power.at.05(sim.prms[,12])

    # Newline for prettier prints
    cat("\n")
}

# Create dataframe for plotting power at various B3s
df.plot <- data.frame(B3.Odds.Seq, B3.Powr.Seq)

# Plot Power Analysis
plotLines(df.plot, 'Moderated Logistic Regression Power', './')

