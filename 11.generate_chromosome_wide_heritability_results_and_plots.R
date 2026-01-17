# generate a dta table with all info
library(data.table)
library(stringr)

# Function to extract values from a single .hsq file
parse_hsq <- function(file) {
  x <- readLines(file)

  # Helper to extract numeric value from a line starting with a pattern
  get_val <- function(pattern) {
    line <- x[grep(pattern, x)]
    if (length(line) == 0) return(NA_real_)
    as.numeric(strsplit(line, "\\s+")[[1]][2])
  }

  # extract SE which is in the 3rd column
  get_se <- function(pattern){
        line <- x[grep(pattern, x)]
        if (length(line) == 0) return(NA_real_)
        parts <- strsplit(line, "\\s+")[[1]]
        as.numeric(parts[3])
}
  data.table(
    file = basename(file),
    chr  = str_extract(basename(file), "\\d+"),  # extract chromosome number
    h2_chr_region = get_val("V\\(G1\\)/Vp"),
    h2_chr_region_se = get_se("V\\(G1\\)/Vp"),
    h2_rest   = get_val("V\\(G2\\)/Vp"),
    h2_rest_se = get_se("V\\(G2\\)/Vp"),
    h2_total  = get_val("Sum of V\\(G\\)/Vp"),
    h2_total_se = get_se("Sum of V\\(G\\)/Vp"),
    V_G1      = get_val("^V\\(G1\\)"),
    V_G2      = get_val("^V\\(G2\\)"),
    V_e       = get_val("^V\\(e\\)"),
    Vp        = get_val("^Vp"),
    LRT       = get_val("^LRT"),
    Pval      = get_val("^Pval")
  )
}
# Read all .hsq files in a directory
hsq_files <- list.files(path = "/hdd/8TBSSD/Steffi/Sheep_genotypes/all_IDs_2583/ssGWAS/Birthweight/gcta_chr_RHM_jan2026/", pattern = "_RHM\\.hsq$", full.names = TRUE)

# Parse all files and combine
results <- rbindlist(lapply(hsq_files, parse_hsq), use.names = TRUE, fill = TRUE)

# Order by chromosome
results[, chr := as.integer(chr)]
results <- results[order(chr)]

# Save to file
fwrite(results, "gcta_chr_RHM_jan2026/RHM_chromosome_wise_all_results.tsv", sep = "\t")

# Plot
library(data.table)
library(ggplot2)

# Load your combined RHM results
rhm <- fread("gcta_chr_RHM_jan2026/RHM_chromosome_wise_all_results.tsv")

# Ensure chromosome is numeric and ordered
rhm[, chr := as.integer(chr)]
setorder(rhm, chr)
# Plot
png("gcta_chr_RHM_jan2026/Birthweight_chr_rhm_heritability.png")
ggplot(rhm, aes(x = chr, y = h2_chr_region)) +
  geom_point(size = 3, colour = "#1f78b4") +
  geom_line(colour = "#1f78b4") +
  theme_minimal(base_size = 14) +
  labs(
    x = "Chromosome",
    y = "Regional heritability (h²)",
    title = "Regional heritability across chromosomes"
  ) +
  scale_x_continuous(breaks = rhm$chr) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

#print(p)
dev.off()

# Highlight significant heritability on chromosomes
png("gcta_chr_RHM_jan2026/Birthweight_chr_rhm_sig_chr_with_heritability.png")
rhm[, sig := Pval < 0.05]
ggplot(rhm, aes(x = chr, y = h2_chr_region, colour = sig)) +
  geom_point(size = 3) +
  geom_line() +
  scale_colour_manual(values = c("FALSE" = "grey60", "TRUE" = "red")) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Chromosome",
    y = "Regional heritability (h²)",
    colour = "Significant",
    title = "Regional heritability across chromosomes (highlighting significant regions)"
  )
dev.off()

