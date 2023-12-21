# Determine percent of variation associated with each PC 

run_PCnum <- function(sce,reduction = "pca") { #, outfile
  pct<-sce[[reduction]]@stdev/sum(sce[[reduction]]@stdev)*100
  # Calculate cumulative percents for each PC
  cumu <- cumsum(pct)
  # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
  co1 <- which(cumu > 90 & pct < 5)[1]
  # Determine the difference between variation of PC and subsequent PC
  # last point where change of % of variation is more than 0.1%.
  co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.05), decreasing = T)[1] + 1
  # Minimum of the two calculation
  pcs <- min(co1, co2)
  # Create a dataframe with values
  # plot_df <- data.frame(pct = pct,cumu = cumu, rank = 1:length(pct))
  # write.csv(plot_df,file = outfile)
  return(pcs)
}
