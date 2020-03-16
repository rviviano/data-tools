# Plot means with 95% confidence intervals for continuous variables and 
# proportions with 95% confidence intervals for categorical variables.
# Grouping variable can be a latent class definition, which is the assumption
# that was made when I wrote these functions.
#
# The functions are designed to be reusable in any script; modifications
# should only need to occur in the main() function
#
# Raymond Viviano
# rayviviano@gmail.com
# October 24th, 2019

# For debugging, uncomment if sourcing script within an R prompt to debug
#browser()

# TODO: Rewrite functions so that you can pass dataframes by reference rather 
#       than by value to improve efficiency. Or switch from dataframes to 
#       data.tables for working with larger datasets.

# TODO: Allow arbitrary data files as input. Right now the script expects a
#       Stata .dta file as a command-line argument.

# TODO: Add easy functionality for changing color palletes 

# TODO: Split up plotting lines with (p <- p + ...) and comments in the 
#       means and proportions functions for better readability/maintenance.

# TODO: Remove wesanderson library as a dependency

# Dependencies
suppressPackageStartupMessages(library(foreign))
suppressPackageStartupMessages(library(survey))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(wesanderson))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(lattice))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(optparse))

# R stack traces are the worst, this option helps error hunting marginally 
options(show.error.locations=TRUE)
options(error=traceback)


#' Convert decimal to percentage
decimal2percent <- function(x, digits = 2, format = "f", ...) {
    paste0(formatC(100 * x, format = format, digits = digits, ...), "%")
}


#' Compute standard error of mean with NA removal support
stderr <- function(x, na.rm=FALSE){
     if(na.rm==TRUE) x <- na.omit(x)
     sd(x)/sqrt(length(x))
}


#' Generates a binary subpopulation variable for a dataframe based on
#' the level of a categorical variable passed to the function.
#' Inputs: 
#'     df:      Dataframe with class var and independent var of interest.
#'     spName:  Binary subpopulation variable to generate
#'     catVar:  Categorical variable to base the subpopulation on 
#'     value:   Value of catVar where spName <- 1
#' Output:
#'     df: Same dataframe with added column
genSubpopFromCategory <- function(df, spName, catVar, value){
    df[, spName] <- 0
    for (i in 1:nrow(df)){
        if (df[i, catVar] == value){
            df[i, spName] <- 1
        }
    }
    return(df)
}


#' Takes a dataframe and a column name of the variable to compute the 
#' Z-scores for. If subpop is specified, only compute z-scores over
#' the specified subpopulation and put NA values in for individuals 
#' that do not belong to the subpopulation. NOTE: R handles complex survey
#' data, namely standard deviations, differently from STATA. Results will not
#' be comparable, use STATA to generate Z Scores if possible.
#'
#' Inputs: 
#'     df:      Complete dataframe 
#'     x:       variable to compute zscore for 
#'     weight:  Participant survey weight
#'     strat:   Stratification
#'     psu:     Primary sampling unit ID
#'     subpop:  (Optional) Subpopulation to compute zscore for, assumes a 
#'              binary variable and computes z-scores based on individuals
#'              with a value of 1.
#' Output:
#'     df:      Same dataframe with added column, columname will be 
#'              <x>_zscore_<subpop> if a subpop was specified,
#'              else, the column name will be <x>_zscore
#'
#' Currently unoptimized for super large dataframes as R passes by value
#' rather than by reference when the dataframe is manipulated within the 
#' function. For [very] large datasets, calling this function a lot could 
#' lead to heap thrashing.
genComplexSurveyZScore <- function(df, x, weight, strat, psu, subpop=NULL){
    # Define survey design
    design <- svydesign(id=as.formula(paste0("~", psu)), 
                        weight=as.formula(paste0("~", weight)), 
                        strata=as.formula(paste0("~", strat)), data=df)

    # Calculate mean and variance accounting for complex survey design
    if (!is.null(subpop)){
        meanSE <- svymean(as.formula(paste0("~", x)), 
                          subset(design, get(subpop)==1), 
                          na.rm=TRUE)

        variance <- svyvar(as.formula(paste0("~", x)), 
                           subset(design, get(subpop)==1), 
                           na.rm=TRUE)
    } else {
        meanSE <- svymean(as.formula(paste0("~", x)), 
                          na.rm=TRUE)

        variance <- svyvar(as.formula(paste0("~", x)), 
                           na.rm=TRUE)
    }

    # Extract mean from vector holding mean & standard error
    meanX <- as.numeric(meanSE[1])

    # Compute standard deviation from variance
    stdevX <- sqrt(as.numeric(variance[1]))

    # TODO: Super inefficient...
    # Create zscore column in the dataframe
    if (!is.null(subpop)){
        df[, paste0(x, "_zscore_", subpop)] <- NA
        for (i in 1:nrow(df)){
            if (df[i, subpop] == 1){
                df[i, paste0(x, "_zscore_", subpop)] <- (df[i, x] - meanX)/stdevX
            }
        }
    } else {
        df[, paste0(x, "_zscore")] <- NA
        for (i in 1:nrow(df)){
            df[i, paste0(x, "_zscore")] <- (df[i, x] - meanX)/stdevX
        }
    }

    # Return dataframe with added zscore column
    return(df)
}


#' Calculate mean and stderr by class
#'
#' Inputs: 
#'     df:      Dataframe with class var and independent var of interest.
#'     cVar:    Class variable specifying the specifc latent class solution 
#'     iVar:    Independent variable
#' Output:
#'     dfStats: Dataframe with mean and stderr for each class
unweightedMeanSE<- function(df, cVar, iVar){
    # Reduce dataframe to only the variables of interest
    df <- df[c(cVar, iVar)]

    # Remove rows with NA
    df <- df[complete.cases(df),]

    # Get number of groups
    numGrps <- max(df[cVar])
    grps  <- NULL
    means <- NULL
    sterr <- NULL

    # Iterate over each group
    for (i in 1:numGrps){
        # Subset the data based on the group
        dfTemp <- df 
        dfTemp <- dfTemp[which(dfTemp[cVar]==i), ]
        # Append group, mean, and standard error to respective vectors
        grps  <- c(grps, i)
        means <- c(means, mean(dfTemp[[iVar]]))
        sterr <- c(sterr, stderr(dfTemp[[iVar]])) 
    }

    # Create a summary dataframe to return to plotting function
    dfStats <- data.frame(Class=as.factor(grps), mean=means, se=sterr)

    return(dfStats)
} 


#' Calculate mean and stderr by class taking into account complex 
#' survey designs. Must specify weight, stratification, and clustering
#'
#' Inputs: 
#'     df:      Dataframe with class var and independent var of interest.
#'     cVar:    Class variable specifying the specifc latent class solution 
#'     iVar:    Independent variable
#'     weight:  Participant survey weight
#'     strat:   Stratification
#'     psu:     Primary sampling unit ID
#' Output:
#'     dfStatsC: Dataframe with mean and stderr for each class (complex)
weightedMeanSE <- function(df, cVar, iVar, weight, strat, psu){
    # Reduce dataframe to only the variables of interest
    df <- df[c(cVar, iVar, weight, strat, psu)]

    # Change column names for easier reference within this function
    colnames(df)[colnames(df)==cVar]   <- "cls"
    colnames(df)[colnames(df)==iVar]   <- "ind"
    colnames(df)[colnames(df)==weight] <- "wgt"
    colnames(df)[colnames(df)==strat]  <- "strt"
    colnames(df)[colnames(df)==psu]    <- "psu"

    # Define survey design
    design <- svydesign(id=~psu, weight=~wgt, strata=~strt, data=df)

    # Calculate mean and standard error
    dfStatsC <- svyby(~ind, ~cls, design, svymean, na.rm=TRUE)
    dfStatsC$Class <- as.factor(dfStatsC$cls)

    return(dfStatsC)

}


#' Calculate proportions of individuals reporting each level of a 
#' categorical variable by class. Also calculate standard errors
#'
#' Inputs: 
#'     df:      Dataframe with class var and independent var of interest.
#'     cVar:    Class var specifying the specific latent class solution 
#'     iVar:    Independent variable
#' Output:
#'     dfStats: Dataframe with proportions of individuals at each level
#'              of the independent variable by class. Includes stderrs.
unweightedProportions <- function(df, cVar, iVar){
    # Reduce dataframe to only the variables of interest
    df <- df[c(cVar, iVar)]

    # Remove rows with NA
    df <- df[complete.cases(df),]

    # Get number of groups
    numGrps <- max(df[cVar])

    # Force categorical variable into factor datatype
    df$iFactor <- as.factor(df[[iVar]])

    # Get levels of independent variable
    lvls <- levels(df$iFactor)

    # Initialize vectors as null type 
    grps   <- NULL
    lvlVec <- NULL
    props  <- NULL
    sterr  <- NULL

    # Iterate over each group
    for(i in 1:numGrps){
        # Subset the data based on the group
        dfSub <- df[which(df[cVar]==i), ]
        # Iterate over each level of the independent variable
        for (j in lvls){
            # Subset dfSub by categorical variable level (subset of subset)
            dfSub2 <- dfSub[which(dfSub$iFactor==j), ]
            # Calculate proportion and accompanying standard error
            proportion <- nrow(dfSub2)/nrow(dfSub)
            standerror <- sqrt(proportion*(1-proportion)/nrow(dfSub))
            # Append group, proportion, and std error to respective vectors
            grps   <- c(grps, i)
            lvlVec <- c(lvlVec, j)
            props  <- c(props, proportion)
            sterr  <- c(sterr, standerror) 
        }
    }

    # Create a summary dataframe to return to plotting function
    dfStats <- data.frame(Class=as.factor(grps), ind=lvls, prop=props, se=sterr)
    
    return(dfStats)
}


#' Calculate proportions of individuals reporting each level of a 
#' categorical variable by class, taking into account complex 
#' survey designs. Also compute stderrs. Must specify weight, 
#' stratification, and clustering
#'
#' Inputs: 
#'     df:      Dataframe with class var and independent var of interest.
#'     cVar:    Class variable specifying the specifc latent class solution 
#'     iVar:    Independent variable
#'     weight:  Participant survey weight
#'     strat:   Stratification
#'     psu:     Primary sampling unit ID
#' Output:
#'     dfStatsC: Dataframe with mean and stderr for each class (complex)
weightedProportions <- function(df, cVar, iVar, weight, strat, psu){
    # Reduce dataframe to only the variables of interest
    df <- df[c(cVar, iVar, weight, strat, psu)]

    # Change column names for easier reference within this function
    colnames(df)[colnames(df)==cVar]   <- "cls"
    colnames(df)[colnames(df)==iVar]   <- "ind"
    colnames(df)[colnames(df)==weight] <- "wgt"
    colnames(df)[colnames(df)==strat]  <- "strt"
    colnames(df)[colnames(df)==psu]    <- "psu"

    # Get number of groups
    numGrps <- max(as.integer(df$cls), na.rm=TRUE)

    # Force categorical variable into factor datatype
    df$ind <- as.factor(df$ind)

    # Define survey design
    des <- svydesign(id=~psu, weight=~wgt, strata=~strt, data=df)

    # Initialize output dataframe
    dfStatsC <- data.frame(Class=integer(), ind=character(), prop=double(), 
                           se=double(), stringsAsFactors=FALSE)

    # Iterate over each group
    for(i in 1:numGrps){
        # Subset the data based on class, get props and SEs
        propsTable <- svymean(~ind, subset(des, cls==i), na.rm=TRUE)
        props <- as.data.frame(propsTable)
        rownames <- row.names(props)
        # Add class, factor level (as character), proportion, and se to df
        for (j in rownames){
            row <- list(i, as.character(substr(j,4,nchar(j))), 
                        props[j,"mean"], props[j,"SE"])
            row <- as.data.frame(row, stringsAsFactors=FALSE)
            colnames(row) <- c("Class", "ind", "prop", "se")
            dfStatsC <- rbind(dfStatsC, row)
        }   
    }
    # Convert columns to factor type
    dfStatsC$Class <- as.factor(dfStatsC$Class)
    dfStatsC$ind <- as.factor(dfStatsC$ind)  

    return(dfStatsC)
}


#' Barplots comparing continuous variables across latent classes
#' Includes means and 95% confidence intervals
#'
#' Inputs: 
#'     df:      Dataframe with class var and independent var of interest.
#'     cVar:    Class colname specifying the specific latent class solution 
#'              to use as the grouping variable for the barplots.
#'     iVar:    Independent variable colname
#'     title:   Title above chart 
#'     outDir:  Directory to save plot to
plotUnweightedMeans <- function(df, cVar, iVar, title, outDir){
    # Calculate mean and standard error for each group
    dfStats <- unweightedMeanSE(df, cVar, iVar)

    # Create barplot of class means with 95% confidence intervals
    p <- ggplot(data=dfStats, aes(x=Class, y=mean, color=Class, group=Class, 
         fill=Class)) + geom_bar(stat='identity', alpha=0.7) + 
         geom_errorbar(aes(ymin = mean - 1.96*se, ymax = mean + 1.96*se), 
         width=.05, alpha=.4, color="black") + geom_text(aes(x=Class, y=mean+2.3*se+mean/20, label=gsub("0\\.", "\\.", round(mean, 2))), 
         color='black', size=4, vjust="inward") + scale_y_continuous(expand=c(0,0)) + 
         theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
         panel.background=element_blank(), axis.line=element_line(colour="black")) +
         labs(y='Mean with 95% Confidence Interval', x='Class') +
         scale_fill_manual(values=wes_palette(name="Moonrise2"), name="Class") + 
         scale_color_manual(values=rep(c("grey"), nlevels(dfStats$Class)), guide=FALSE) +
         ggtitle(paste("Mean ", title, " Unweighted", sep=""))

    # Save Plots
    pngOut <- paste(outDir, "/", cVar, "-", iVar, ".png", sep="")    
    svgOut <- paste(outDir, "/", cVar, "-", iVar, ".svg", sep="")  
    
    ggsave(pngOut, p)
    ggsave(svgOut, p)
}


#' Barplots comparing continuous variables across latent classes
#' Includes means and 95% confidence intervals
#'
#' Inputs: 
#'     df:      Dataframe with class var and independent var of interest.
#'     cVar:    Class colname specifying the specific latent class solution 
#'              to use as the grouping variable for the barplots.
#'     iVar:    Independent variable colname
#'     title:   Title above chart 
#'     wgt:     Participant survey weight
#'     strat:   Stratification
#'     psu:     Primary sampling unit ID
#'     outDir:  Directory to save plot to
plotWeightedMeans <- function(df, cVar, iVar, wgt, strat, psu, title, outDir){
    # Calculate mean and standard error for each group (complex)
    dfStatsC <- weightedMeanSE(df, cVar, iVar, wgt, strat, psu)

    # Create barplot of class means with 95% confidence intervals
    p <- ggplot(data=dfStatsC, aes(x=Class, y=ind, color=Class, group=Class, 
         fill=Class)) + geom_bar(stat='identity',alpha=0.7) + 
         geom_errorbar(aes(ymin = ind - 1.96*se, ymax = ind + 1.96*se), 
         width=.05, alpha=.4, color="black") + geom_text(aes(x=Class, y=ind+2.3*se+ind/20, label=gsub("0\\.", "\\.", round(ind, 2))), 
         color='black', size=4, vjust="inward") + scale_y_continuous(expand=c(0,0)) + 
         theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
         panel.background=element_blank(), axis.line=element_line(colour="black")) +
         labs(y='Mean with 95% Confidence Interval', x='Class') +
         scale_fill_manual(values=wes_palette(name="Moonrise2"), name="Class") + 
         scale_color_manual(values=rep(c("grey"), nlevels(dfStatsC$Class)), guide=FALSE) + 
         ggtitle(paste("Mean ", title, " Accounting for Complex Survey Design", sep=""))

    # Save Plots
    pngOut <- paste(outDir, "/", cVar, "-", iVar, "-cmplxsrvy.png", sep="")    
    svgOut <- paste(outDir, "/", cVar, "-", iVar, "-cmplxsrvy.svg", sep="")  
    
    ggsave(pngOut, p)
    ggsave(svgOut, p)
}


#' Barplots comparing categorical variable proportions across latent classes
#' Includes proportions and 95% confidence intervals
#'
#' Inputs: 
#'     df:      Dataframe with class var and independent var of interest.
#'     cVar:    Class variable specifying the specifc latent class solution 
#'              to use as the grouping variable for the barplots.
#'     iVar:    Independent variable
#'     title:   Title above chart 
#'     outDir:  Directory to save plot to
plotProportions <- function(df, cVar, iVar, title, outDir){
    # Calculate proportion and standard error for each group
    dfStats <- unweightedProportions(df, cVar, iVar)

    # Create barplot of proportions with 95% confidence intervals
    p <- ggplot(data=dfStats, aes(x=Class, y=prop, color=ind, 
         fill=ind)) + geom_bar(stat='identity', alpha=0.85, width=.79, position=position_dodge(.8)) + 
         geom_errorbar(aes(ymin = prop - 1.96*se, ymax = prop + 1.96*se), 
         width=.05, alpha=0.4, color="black", position = position_dodge(0.8)) + 
         geom_text(aes(x=Class, y=dfStats$prop+2*se+.015, label=gsub("0\\.", "\\.", round(prop, 2))), 
         color='black', position = position_dodge(0.8), size=2.7, vjust="inward") + scale_y_continuous(expand=c(0,0)) + 
         theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
         panel.background=element_blank(), axis.line=element_line(colour="black"), legend.title=element_blank()) +
         labs(y='Proportion with 95% Confidence Interval', x='Class') +
         scale_fill_manual(values=brewer.pal(n=nlevels(dfStats$ind), name="BrBG"), name=title) + 
         scale_color_manual(values=rep(c("grey"), nlevels(dfStats$ind)), guide=FALSE) + 
         ggtitle(paste(title, " Unweighted Proportions", sep=""))

    # Save Plots
    pngOut <- paste(outDir, "/", cVar, "-", iVar, ".png", sep="")    
    svgOut <- paste(outDir, "/", cVar, "-", iVar, ".svg", sep="")  
    
    ggsave(pngOut, p)
    ggsave(svgOut, p)
}


#' Barplots comparing categorical variable proportions across latent classes
#' Includes proportions and 95% confidence intervals
#'
#' Inputs: 
#'     df:      Dataframe with class var and independent var of interest.
#'     cVar:    Class variable specifying the specifc latent class solution 
#'              to use as the grouping variable for the barplots.
#'     iVar:    Independent variable
#'     title:   Title above chart 
#'     wgt:     Participant survey weight
#'     strat:   Stratification
#'     psu:     Primary sampling unit ID
#'     outDir:  Directory to save plot to
plotWeightedProps <- function(df, cVar, iVar, wgt, strat, psu, title, outDir){
    # Calculate proportion and standard error for each group
    dfStatsC <- weightedProportions(df, cVar, iVar, wgt, strat, psu)

    # Create barplot of proportions with 95% confidence intervals
    p <- ggplot(data=dfStatsC, aes(x=Class, y=prop, color=ind, 
         fill=ind)) + geom_bar(stat='identity', alpha=0.85, width=.79, position=position_dodge(.8)) + 
         geom_errorbar(aes(ymin = prop - 1.96*se, ymax = prop + 1.96*se), 
         width=.05, alpha=0.4, color="black", position = position_dodge(0.8)) + 
         geom_text(aes(x=Class, y=dfStatsC$prop+2*se+.015, label=gsub("0\\.", "\\.", round(prop, 2))), 
         color='black', position = position_dodge(0.8), size=2.7, vjust="inward") + scale_y_continuous(expand=c(0,0)) + 
         theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), 
         panel.background=element_blank(), axis.line=element_line(colour="black"), legend.title=element_blank()) +
         labs(y='Proportion with 95% Confidence Interval', x='Class') +
         scale_fill_manual(values=brewer.pal(n=nlevels(dfStatsC$ind), name="BrBG"), name=title) + 
         scale_color_manual(values=rep(c("grey"), nlevels(dfStatsC$ind)), guide=FALSE) + 
         ggtitle(paste(title, " Proportions Accounting for Complex Survey Design", sep=""))

    # Save Plots
    pngOut <- paste(outDir, "/", cVar, "-", iVar, "-cmplxsrvy.png", sep="")    
    svgOut <- paste(outDir, "/", cVar, "-", iVar, "-cmplxsrvy.svg", sep="")  
    
    ggsave(pngOut, p)
    ggsave(svgOut, p)
}


#' Lineplot comparing continuous variable zscore means across latent classes
#' Includes proportions and 95% confidence intervals. This information
#' could be conveyed with a barplot, but it's much cleaner to group
#' multiple variables together with this form of zscore barplot
#'
#' Inputs: 
#'     df:       Dataframe with class var and independent vars of interest.
#'     cVar:     Class variable specifying the specifc latent class solution 
#'               to use as the grouping variable for the line plot.
#'     iList:    List of independent variables to plot
#'     axisLbls: Labels that correspond to iList variables for x-axis
#'     wgt:      Participant survey weight
#'     strat:    Stratification
#'     psu:      Primary sampling unit ID
#'     title:    Title above chart 
#'     clrVals:  Vector of hex vals for coloring, e.g. c("#999999",... etc.) 
#'     outDir:   Directory to save plot to
plotZScoreLinePlot <- function(df, cVar, iList, axisLbls, wgt, strat,
                               psu, title, clrVals, outDir){
    # Initialize vectors as null type to hold all means/stderrors
    grps   <- NULL
    lblVec <- NULL
    means  <- NULL
    sterr  <- NULL

    # Calculate mean and standard error for each variable over each group
    for (i in 1:length(iList)){
        iVar <- iList[i]
        iLbl <- axisLbls[i]
        # Get mean and standard error of z-score for each class 
        dfStatsC <- weightedMeanSE(df, cVar, iVar, wgt, strat, psu)
        # Loop over the classes
        for (j in levels(dfStatsC$Class)){
            # Append class, axis lbl, mean, and stderror to respective vectors
            grps   <- c(grps, j)
            lblVec <- c(lblVec, iLbl)
            means  <- c(means, dfStatsC[j, "ind"])
            sterr  <- c(sterr, dfStatsC[j, "se"])
        }
    }

    # Create a dataframe specific to plotting
    dfPlot <- data.frame(Class=as.factor(grps), lbls=lblVec, 
                         means=means, se=sterr, stringsAsFactors=FALSE)

    # Create plot of means with 95% confidence intervals - Main Layer
    p <- ggplot(data=dfPlot, aes(x=lbls, y=means, color=Class, 
                group=Class, fill=Class)) 

    # Add lines for each group
    p <- p + geom_line(aes(linetype=Class), stat='identity', 
                       size=.85, alpha=0.9) 
    
    # Unique shape for each group placed at each data point
    p <- p + geom_point(aes(shape=Class), size=2.5, alpha=.45) 

    # Add 95% Confidence Intervals
    p <- p + geom_errorbar(aes(ymin = means - 1.96*se, ymax = means + 1.96*se), 
                           width=.05, alpha=.6, color="black") 


    # Pad the distance between the x-axis and smallest plotted value slightly
    p <- p + scale_y_continuous(expand=c(0.05, 0))

    # Pad/add extra space to the ends of the x-axis
    p <- p + scale_x_discrete(expand=c(0.05, 0.05)) 

    # Remove gray background, format x-axis text, format legend 
    p <- p + theme(panel.grid.major=element_blank(), 
                   panel.grid.minor=element_blank(), 
                   panel.background=element_blank(), 
                   axis.line=element_line(colour="black"), 
                   axis.text.x=element_text(angle=50, size=8, hjust=1, vjust=1), 
                   legend.key=element_rect(color=NA, fill=NA),
                   legend.background=element_blank(),
                   legend.key.width=unit(3, "line"),
                   legend.justification=c(1, 1),
                   legend.position=c(1, 1)) 
    
    # Add y-axis label
    p <- p + labs(y='Mean Z-Scores with 95% CI', x=NULL) 

    # Set the color scheme for the plot
    p <- p + scale_fill_manual(values=clrVals, name="Class", guide=FALSE) 
    p <- p + scale_color_manual(values=clrVals, name="Class") 

    # Remove any components to the legend that could relate to alpha or size
    p <- p + scale_alpha(guide='none') 
    p <- p + scale_size(guide='none') 

    # Add a title to the plot
    p <- p + ggtitle(paste0("Mean ", title, " Accounting for Survey Design")) 

    # Set the x/y axis ratio based on number of independent variables
    p <- p + coord_fixed(ratio=12/length(axisLbls))

    # Save Plots
    pngOut <- paste(outDir, "/", cVar, "-zscoreline-cmplxsrvy.png", sep="")    
    svgOut <- paste(outDir, "/", cVar, "-zscoreline-cmplxsrvy.svg", sep="")  
    
    ggsave(pngOut, p)
    ggsave(svgOut, p)
}


#' Calls to all plotting functions pulled out of the main function to 
#' "hide" from user modification. All user modifications based on the 
#' project should occur in the main function.
plotAllBar <- function(df, classVars, varsContinuous, titlesContinuous, 
                       varsCategorical, titlesCategorical, figureDir){
    # Plot continuous variables over class solutions (unweighted)
    for (i in 1:length(classVars)){
        for(j in 1:length(varsContinuous)){
            cVar <- classVars[i]
            iVar <- varsContinuous[j]
            title <- titlesContinuous[j]
            print(toString(c("Plotting Mean", title, "Unweighted", cVar)))
            plotUnweightedMeans(df, cVar, iVar, title, figureDir)
        }
    }

    # Plot continuous variables over class solutions (complex survey)
    for (i in 1:length(classVars)){
        for(j in 1:length(varsContinuous)){
            cVar <- classVars[i]
            iVar <- varsContinuous[j]
            title <- titlesContinuous[j]
            print(toString(c("Plotting Mean", title, "Accounting for Survey Design", cVar)))
            plotWeightedMeans(df, cVar, iVar, "weight_final_norm_overall", 
                              "strat", "psu_id", title, figureDir)
        }
    }

    # Plot categorical variables over class solutions (unweighted)
    for (i in 1:length(classVars)){
        for(j in 1:length(varsCategorical)){
            cVar <- classVars[i]
            iVar <- varsCategorical[j]
            title <- titlesCategorical[j]
            print(toString(c("Plotting", title, "Unweighted Proportions", cVar)))
            plotProportions(df, cVar, iVar, title, figureDir)
        }
    }

    # Plot categorical variables over class solutions (complex survey)
    for (i in 1:length(classVars)){
        for(j in 1:length(varsCategorical)){
            cVar <- classVars[i]
            iVar <- varsCategorical[j]
            title <- titlesCategorical[j]
            print(toString(c("Plotting", title, "Proportions Accounting for Survey Design", cVar)))
            plotWeightedProps(df, cVar, iVar, "weight_final_norm_overall", 
                              "strat", "psu_id", title, figureDir)
        }
    }
}


#' setOptionList
#' Create viable options for running this script
setOptionList <- function(){
    optionList = list(make_option(c("-i", "--input"), type="character",
                                  default=NULL, help=paste0("File input (",
                                  ".dta of main dataset).")),
                      make_option(c("-o", "--output"), type="character",
                                  default=NULL, help="Output figure directory"))
    return(optionList)
}


#' verifyArgs
#' Takes a list of options and verifies that they are correctly specified
#' Also takes the optParser object to print help msgs if needed.
#' Stops the script if an argument is invalid
verifyArgs <- function(opt, optParser){ 
    # Check that a file is input at all
    if (is.null(opt$input)){
        print_help(optParser)
        stop("Must enter an input file (.dta of main dataset with LPA solutions).")
    }

    # Check the existence of the input file path
    if (!(file.exists(opt$input))){
        print_help(optParser)
        stop("Input file does not exist. Please check file path.")
    }

    # Check that an output argument is specified
    if (is.null(opt$output)){
        print_help(optParser)
        stop("Must enter a figure ouput directory.")
    }

    # Check that the ouput directory exists
    if (!dir.exists(opt$output)){
        print_help(optParser)
        stop("Figure output directory must exist. Check filepath.")
    }
}


#' Main script logic, all changes for future projects should occur within
#' this function only. Specify paths to your dataset with 'datapath' and 
#' to the directory that you would like to output figures to with 
#' 'figureDir.' The script currently assumes that the data is in a Stata 
#' .dta. If not loading from a stata .dta, change the 'read.dta' line to fit
#' your data. Simply specify the grouping or class variables, the continuous 
#' variables, and the categorical variables. Also create separate vectors of 
#' titles for the continuous and categorical variables. The script should 
#' be able to create reasonable bar plots for means and propotions that
#' can be easily modified for other projects
main <- function(){
    # Set up the options parser
    optionList <- setOptionList()
    optParser = OptionParser(option_list=optionList);

    # Get command line arguments
    opt = parse_args(optParser);

    # Verify arguments
    verifyArgs(opt, optParser)

    # Path to main dataset with LPA solutions already merged in
    dataPath <- opt$input

    # Load dataset for plotting
    df <- read.dta(dataPath)   

    # Location to save plots to
    figureDir <- opt$output

    # NOTE: R Survey computes different standard deviation compared to STATA
    #       do not use the genComplexSurveyZScore over STATA.
 
}


### Run script ###
main()