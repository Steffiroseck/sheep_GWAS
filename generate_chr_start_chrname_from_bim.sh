awk '$1==26 {print $1"\t"$4"\t"$2}' 2583_UKIDs_afterQC_pruned_Hapmap_imputed_pruned.bim | sort -k2,2n > chr26.map
