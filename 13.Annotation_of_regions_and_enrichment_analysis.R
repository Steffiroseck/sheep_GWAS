# Load the necessary libraries
library(qqman)
library(stringr)
library(dplyr)
library(ggplot2)
library(purrr)
library(data.table)

# Define the base folder
base_dir <- "age_at_slaughter_154/RHM_results/"

# List all files (non-recursive) in the folder
all_files <- list.files(base_dir, full.names = TRUE)

# Filter files that contain the desired patterns and have the correct extension
target_files <- all_files[grepl("Age_at_slaughter(_RHM)?", all_files) & grepl("\\.Results$", all_files)]
target_files


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

clean_results_files("age_at_slaughter_154/RHM_results/")

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
  png(file.path(base_dir, paste0(main_topic, "_manhattan.png")), width = 800, height = 600)
  manhattan(data, chr = "Chr", bp = "End", snp = "Start", p = "Pval", ylim = c(0, 8),
            main = main_topic, cex = 1.5, cex.axis = 0.9,
            col = c("blue", "red"), suggestiveline = suggest, genomewideline = genomewide,
            chrlabs = as.character(c(1:22, 24:26)))
  dev.off()  # Close the png device
  
  # Calculate lambda for genomic control
  lambda1 <- qchisq(1 - median(data$Pval), 1) / qchisq(0.5, 1)
  
  # QQ plot
  observed <- -log10(sort(data$Pval, decreasing = FALSE))
  expected <- -log10(ppoints(length(observed)))
  lambda2 <- median(observed) / median(expected) # λGC > 1.05 suggests stratification
  
  png(file.path(base_dir, paste0(main_topic, "_qq.png")), width = 800, height = 600)
  qq(data$Pval, main = main_topic)  # Use cleaned data
  mtext(paste("GIF =", round(lambda2, 2)), side = 3, line = 0.5, adj = 1)
  dev.off()  # Close the png device
}


\ Extract the region with significant association found
rhm_10_5 <- read.csv("RHM_results/Age_at_slaughter_RHM_10_5.Results", sep="\t")
main_topic <- "Age at slaughter RHM_10_5"
# Total SNPs after QC
total_snps <- 43823

# Window size and overlap
window_size <- 10
overlap <- 5

# Estimate number of windows
num_windows <-  floor((total_snps - window_size) / overlap) + 1
num_windows

# Calculate genome-wide significance threshold
genomewide <- -log10(0.05 / num_windows)
genomewide
suggest <- -log10(1 / num_windows)
suggest

manhattan(rhm_10_5, chr = "Chr", bp = "End", snp = "Start", p = "Pval", ylim = c(0, 6),
            main = main_topic, cex = 0.9, cex.axis = 0.9,
            col = c("blue", "red"), suggestiveline = suggest, genomewideline = genomewide,
            chrlabs = as.character(c(1:22, 24:26)))


# For finding significant regions associated
genomewide <- 0.05 / num_windows
suggest <- 1 / num_windows
sig_genome_wide <- rhm_10_5[rhm_10_5$Pval < genomewide, ]
sig_suggestive  <- rhm_10_5[rhm_10_5$Pval < suggest, ]
suggest

# Extract genes

# ###############################################
# # For single mlma file
# ################################################
#
all_gwas <-rhm_10_5
# head(all_gwas)


# filter significant SNPs based on p value, here we use first row with lowest pvalue
#sig_snps <- all_gwas %>% filter(P.adj < 0.4)
all_gwas_sort <- all_gwas[order(all_gwas$Pval, decreasing = FALSE), ]
sig_snps <- all_gwas_sort[1, , drop=FALSE]
dim(sig_snps)
sig_snps

#write.csv(sig_snps, "Birth_weight/gcta_results/BIRTH_WEIGHT_CLEAN_INT_significant_SNPs.csv", row.names=FALSE)

# map SNPs to genes
library(biomaRt)
ensembl <- useMart("ensembl", dataset = "oaries_gene_ensembl")

# map the start and end snp from the SNP.list26_71_80 file of age at slaughter. Then start position will be the sirst snp in that file's position and end position will be last snp in that file. make sure to consider that overlaps will be present as we analysed 10 SNP with 5 snp overlap. You can use SNP_annotation file to get snp positions.
start_position <- 5333971 
end_position <- 5770772
regions <- paste0(as.character(sig_snps$Chr), ":", start_position - 50000, ":", end_position + 50000)
regions

res_all <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "chromosome_name",
                 "start_position", "end_position"),
  filters = "chromosomal_region",
  values = regions,
  mart = ensembl
)
# distance in base pairs
res_all$distance_bp <- 5770772- 5333971
# convert to kilobases
res_all$distance_kb <- res_all$distance_bp / 1000
res_all

write.csv(res_all, "age_at_slaughter_154/Age_at_slaughter_Chromosome26_RHM_genes.csv")
