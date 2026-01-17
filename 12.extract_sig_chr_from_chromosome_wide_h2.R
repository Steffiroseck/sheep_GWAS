setwd("C:/Users/afbi-roses/SHEEP_GENOTYPES/Data_analysed_for_gwas_paper/RHM_all_traits_chromosomewise/")
bw <- read.csv("Birthweight_RHM_chromosome_wise_all_results.tsv", sep="\t")
head(bw)
genomewide <- 0.05 / 26
suggest <- 1 / 26
sig_genome_wide <- bw[bw$Pval < genomewide, ]
sig_suggestive  <- bw[bw$Pval < suggest, ]
