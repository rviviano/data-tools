library(ggplot2)

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
    p <- ggplot(data=df, aes(x=scd_slp_to_obj_slp, y=prop, color=Path, 
                group=Path, fill=Path)) 

    # Add lines for each group
    p <- p + geom_line(aes(linetype=Path), stat='identity', size=.85, alpha=0.9) 
    
    # Unique shape for each group placed at each data point
    p <- p + geom_point(aes(shape=Path), size=2.5, alpha=.45) 

    # Pad the distance between the x-axis and smallest plotted value slightly
    p <- p + scale_y_continuous(limits=c(0,1))

    # Pad/add extra space to the ends of the x-axis
    p <- p + scale_x_continuous(expand=c(0.05, 0.05)) 

    # Reverse x-axis when dealing with negative numbers
    p <- p + scale_x_reverse()

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
                   legend.position=c(1, 1)) 
    
    # Add axis labels
    p <- p + labs(y='Proportion of Simulations with RMSEA<.05 and Path of Interest P<.05', 
                  x='Path Coefficient for Objective Performance Slope ON SCD Slope')  

    # Remove any components to the legend that could relate to alpha or size
    p <- p + scale_alpha(guide='none') 
    p <- p + scale_size(guide='none') 

    # Add a title to the plot
    p <- p + ggtitle(title) 

    # Save Plots
    pngOut <- paste(outDir, "/parallelGrowthPower-linePlot.png", sep="")    
    svgOut <- paste(outDir, "/parallelGrowthPower-linePlot.svg", sep="")   
    
    ggsave(pngOut, p)
    ggsave(svgOut, p)
}


# Path to your csv file
csv_path <- ""

df <- read.csv(csv_path)

temp1 <- df[,c('scd_slp_to_obj_slp', 'rmsea05_si_path05_prop')]
temp2 <- df[,c('scd_slp_to_obj_slp', 'rmsea05_ss_path05_prop')]


colnames(temp1) <- c('scd_slp_to_obj_slp','prop')
temp1[,'Path'] <- "Objective Performance Slope ON SCD Intercept"

colnames(temp2) <- c('scd_slp_to_obj_slp','prop')
temp2[,'Path'] <- "Objective Performance Slope ON SCD Slope"

df <- rbind(temp1, temp2)
df$Path <- as.factor(df$Path)


# # Save data as csv
# write.csv(df, "./parallelGrowth.csv")

plotLines(df, 'Parallel Growth Power Analysis. 100 Control, 100 MCI', './')