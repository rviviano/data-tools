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
    p <- ggplot(data=df, aes(x=s_mci, y=r05p05, color=svar, 
                group=svar, fill=svar)) 

    # Add lines for each group
    p <- p + geom_line(aes(linetype=svar), stat='identity', size=.85, alpha=0.9) 
    
    # Unique shape for each group placed at each data point
    p <- p + geom_point(aes(shape=svar), size=2.5, alpha=.45) 

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
    p <- p + labs(y='Proportion of Simulations with RMSEA<.05 and S~MCI P<.05', 
                  x='S~MCI Unstandardized Path Coefficient')  

    # Remove any components to the legend that could relate to alpha or size
    p <- p + scale_alpha(guide='none') 
    p <- p + scale_size(guide='none') 

    # Add a title to the plot
    p <- p + ggtitle(title) 

    # # Set the x/y axis ratio based on number of independent variables
    # p <- p + coord_fixed(ratio=.4)

    # Save Plots
    pngOut <- paste(outDir, "/conditonalLGMPowerAnalysis-linePlot.png", sep="")    
    svgOut <- paste(outDir, "/conditonalLGMPowerAnalysis-linePlot.svg", sep="")   
    
    ggsave(pngOut, p)
    ggsave(svgOut, p)
}


# Change file load information here
load('./sim_df_m1.Rda')
df <- sim_df_m1
df$svar <- factor(df$svar)

# # Save data as csv
# write.csv(df, "./sim_df_m1.csv")

plotLines(df, 'Conditional LGM Power Analysis. 100 Control, 100 MCI', './')