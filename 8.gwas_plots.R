library(dplyr)
library(ggplot2)
library(qqman)

setwd("C:/Users/afbi-roses/Documents/QUB_PhD/SHEEP_GWAS/")
mlma_files <- list.files(Age_at_slaughter/gcta_results/", pattern = "*.mlma", full.names=TRUE)
process_mlma_file <- function(file_path) {
  # Read in the MLMA file (assuming columns like 'SNP', 'P', and 'Trait')
  gwasdata <- read.table(file_path, header = TRUE)
  
  # Assuming the MLMA file has 'SNP', 'P' (p-value), and 'Trait' column
  colnames(gwasdata)=c("CHR", "SNP", "BP", "A1", "A2", "FREQ", "B", "SE", "P")
  trait_name <- basename(file_path) # Get the trait name from the filename (or extract from column)
  
  print(nrow(gwasdata))
  
  # Genomic Inflation factor
  lambda1=qchisq(1-median(gwasdata$P),1)/qchisq(0.5,1)

  # Genomic control
  observed <- -log10(sort(gwasdata$P, decreasing=FALSE))
  expected <- -log10(ppoints(length(observed)))
  lambda2 <- median(observed) / median(expected)# λGC > 1.05 suggests stratification
  
  # Generate QQ plot
  #qq_plot_file <- paste0("qqplot_", gsub("\\.mlma$", "", trait_name), ".png")
  #png(qq_plot_file)
  qq(gwasdata$P)
  title(main = paste("QQ Plot for", trait_name))
  mtext(paste("GIF =", round(lambda1, 3)), side=3, line = 0.5, adj=1)
  #mtext(paste("GC =", round(lambda2, 3)), side=1, line = 0.5, adj=1)
  #dev.off()

  # Generate Manhattan plot
  #manhattan_plot_file <- paste0("manhattan_", gsub("\\.mlma$", "", trait_name), ".png")
  #png(manhattan_plot_file)
  ## Set threshold for Bonferroni correction
  manhattan(gwasdata, main = paste("Manhattan Plot for", trait_name), ylim = c(0, 6),
            col = c("blue4", "orange3"),
            suggestiveline = -log10(1/nrow(gwasdata)), # Suggestive threshold
            genomewideline = -log10(0.05/nrow(gwasdata))) # Genome-wide significance
            #annotatePval = 0.05) #Annotate top hits
  #dev.off()
  
  print(paste("Generated plots for trait:", trait_name))

  #qqplot for the p-values
  qqnorm(gwasdata$P, pch = 1, frame = FALSE)
  qqline(gwasdata$P, col = "steelblue", lwd = 2)
  points(seq(-4,3.9999,8/length(gwasdata$Bonferroni)) ,sort(gwasdata$Bonferroni), col = "red", cex=0.4,lwd = .2)

}
