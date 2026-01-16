# Load the necessary libraries
library(qqman)
library(stringr)

setwd("C:/Users/afbi-roses/Documents/QUB_PhD/SHEEP_GWAS/")

# Define the base folder
base_dir <- "RHM_results/"

# List all files (non-recursive) in the folder
all_files <- list.files(base_dir, full.names = TRUE)

# Filter files that contain the desired patterns and have the correct extension
target_files <- all_files[grepl("Age_at_slaughter(_RHM)?", all_files) & grepl("\\.Results$", all_files)]
#target_files


# clean the files as some files have less than 13 columns.
clean_results_files <- function(directory) {
  # Get all .results files in the specified directory
  files <- list.files(
    path = directory, 
    pattern = "\\.Results$", 
    full.names = TRUE
  )
  
  for (file in files) {
    # Read all lines from the file
    lines <- readLines(file, warn = FALSE)
    
    # Remove empty lines
    lines <- lines[lines != ""]
    if (length(lines) < 0) next # Skip if file is empty
    
    
    # Identify comment lines (start with '#')
    is_comment <- grepl("^#", lines)
    
    # Find header (first non-comment line)
    header_idx <- which(!is_comment)[1]
    
    # Split lines into columns by whitespace
    cols_list <- strsplit(lines, "\\s+")
    n_cols <- sapply(cols_list, length)
    
    # Initialize vector to mark lines to keep
    keep <- logical(length(lines))
    
    # Always keep comments and header
    keep[is_comment] <- TRUE
    if (!is.na(header_idx)) keep[header_idx] <- TRUE
    
    # Validate data rows (non-comment/non-header): Require exactly 13 columns
    data_rows <- !is_comment & seq_along(lines) != header_idx
    keep[data_rows] <- n_cols[data_rows] == 13
    
    # Write cleaned lines back to the file
    writeLines(lines[keep], file)
  }
}

clean_results_files("RHM_results/")

### For loop
# Loop over each results file
for (file in target_files) {
  # Read the data from the results file
  data <- read.table(file, header = TRUE)
  
  # Extract the filename without the path and extension for naming
  file_name <- basename(file) # Get the base filename
  main_topic <- gsub("\\.Results$", "", file_name) # Remove the .Results extension for the title
  
  # Calculate genome-wide significance threshold
  genomewide <- -log10(0.05 / nrow(data))
  suggest <- -log10(1 / nrow(data))
  
  # Create a Manhattan plot
  #png(file.path(results_dir, paste0(main_topic, "_manhattan.png")), width = 800, height = 600)
  manhattan(data, chr = "Chr", bp = "End", snp = "Start", p = "Pval", ylim = c(0, 6),
            main = main_topic, cex = 0.9, cex.axis = 0.9,
            col = c("blue", "red"), suggestiveline = suggest, genomewideline = genomewide,
            chrlabs = as.character(c(1:22, 24:26)))
  #dev.off()  # Close the png device
  
  # Calculate lambda for genomic control
  lambda1 <- qchisq(1 - median(data$Pval), 1) / qchisq(0.5, 1)
  
  # QQ plot
  observed <- -log10(sort(data$Pval, decreasing = FALSE))
  expected <- -log10(ppoints(length(observed)))
  lambda2 <- median(observed) / median(expected) # λGC > 1.05 suggests stratification
  
  #png(file.path(results_dir, paste0(main_topic, "_qq.png")), width = 800, height = 600)
  qq(data$Pval, main = main_topic)  # Use cleaned data
  mtext(paste("GIF =", round(lambda2, 2)), side = 3, line = 0.5, adj = 1)
  #dev.off()  # Close the png device
}
