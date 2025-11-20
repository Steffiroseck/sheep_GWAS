# Assume that you have run earlier codes and you have variables storing files for this analysis

#In the terminal:
WORKDIR="/hdd/8TBSSD/Steffi/Sheep_genotypes/all_IDs_2583/ssGWAS/Age_at_slaughter"
cd $WORKDIR/
phenotype="$WORKDIR/phenotypes.txt"
sex="$WORKDIR/sex.txt"
COVAR_FILE="covariates.txt"
IMPUTEDIR="/hdd/8TBSSD/Steffi/Sheep_genotypes/"

qc_out="lambs_afterQC_pruned"

OUT_DIR="gcta_results"
mkdir -p "$OUT_DIR"

# generate vcf file from the qc cleaned file
plink --bfile QC/$qc_out \
	--recode vcf \
	--out QC/$qc_out \
	--sheep
	
# Impute genotypes
beagle ref=$IMPUTEDIR/Imputation/SNP50_Breedv1_phased.vcf.gz gt=QC/$qc_out.vcf out=QC/$qc_out"_Hapmap_imputed"

gunzip QC/$qc_out"_Hapmap_imputed".vcf.gz

plink  --vcf QC/$qc_out"_Hapmap_imputed".vcf	--recode --out QC/$qc_out"_Hapmap_imputed" --sheep --allow-extra-chr	--allow-no-sex
plink --file QC/$qc_out"_Hapmap_imputed" --make-bed --out QC/$qc_out"_Hapmap_imputed"	-sheep --allow-extra-chr

# checking QC again
plink --bfile QC/$qc_out"_Hapmap_imputed" \
	--maf 0.05 \
	--geno 0.1 \
	--mind 0.1 \
	--indep-pairwise 50 5 0.8 \
	--allow-extra-chr \
	--allow-no-sex \
	--make-bed \
	--sheep \
	--out QC/$qc_out"_Hapmap_imputed_cleaned" 


# If running pruning in the previous step, then extract the pruned SNPs
plink --bfile QC/$qc_out"_Hapmap_imputed_cleaned" --extract QC/$qc_out"_Hapmap_imputed_cleaned.prune.in" --out QC/$qc_out"_Hapmap_imputed_pruned" --make-bed --sheep --allow-extra-chr


# Update sex info again to the imputed file

for base in *_Hapmap_imputed_pruned.bed; do
  name="${base%.bed}"
  plink --bfile "$name" \
        --update-sex $sex \
        --make-bed \
        --out "${name}" \
        --sheep \
        --allow-extra-chr \
        --allow-no-sex
done

# lets try to add the phenotypes to the imputed.fam file.
nrMeasures=$(head -n 1 "$phenotype" | awk -F' ' '{for(i=3;i<=13;i++) printf $i (i<13 ? " " : "")}')
#Trim extra white spaces;
nrMeasures=$(echo "$nrMeasures" | tr -d '\r' | xargs)
echo "Selected columns:$nrMeasures"
#Now add the phenotypes to .fam file. We will use the phenotypes stored in nrMeaures variable.
for col in $nrMeasures; 
do 
	echo "Updating phenotype for: $col";
	plink --bfile QC/$qc_out"_Hapmap_imputed_pruned" --pheno $phenotype --pheno-name $col --recode --make-bed --out QC/$qc_out"_phenoAdded" --sheep --allow-no-sex --allow-extra-chr; 
	
	echo "PLINK command for $col completed."; 
done

# create the GRM
# Use autosome chromosomes only to avoid sex chromosome bias
gcta64 --bfile QC/$qc_out"_phenoAdded"	--autosome-num 26 --maf 0.05 --make-grm --out QC/$qc_out"_GRM_autosomes"	--thread-num 10


# GCTA GWAS
GRM_PREFIX="QC/$qc_out"_GRM_autosomes""

# Get the trait names from the phenotype file
TRAIT_NAMES=($(head -n 1 "$phenotype" | cut -f3-))
NUM_TRAITS=${#TRAIT_NAMES[@]}

# Loop over each trait
for ((i=1; i<=NUM_TRAITS; i++)); do
    trait=${TRAIT_NAMES[$((i-1))]}
	
	 # Remove carriage return characters if present
     trait=${trait//$'\r'/}
	
    echo "Running GWAS for trait: $trait"
		# Assuming .bed files are named as "<trait_name>.bed"
		BED_FILE="QC/$qc_out"_phenoAdded""

		# Run GCTA MLMA
		gcta64 --grm "$GRM_PREFIX" \
			   --autosome-num 26 \
			   --autosome \
			   --pheno "$phenotype" \
			   --mpheno $i \
			   --mlma \
			   --bfile "$BED_FILE" \
			   --out gcta_results/${trait}"_hapmap"  

		# Create a summary with top SNPs (e.g., P < 1e-5)
		#awk 'NR==1 || $9 < 0.05' "$OUT_DIR/${trait}_gwas.mlma" > "$OUT_DIR/${trait}_top_hits.txt"

done

