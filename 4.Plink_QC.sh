#The terminal () has all the genotype files in .bed, .bim and .fam format. For each sample, these three files will be present. 
#We need to merge all samples into one file for further analysis. To do that, In terminal, 
#we create an output file to store the file names and store it into a variable.

#In the terminal:
WORKDIR="/hdd/8TBSSD/Steffi/Sheep_genotypes/all_IDs_2583/ssGWAS/Age_at_slaughter"
cd $WORKDIR/
output_file="$WORKDIR/Filenames_for_plink_merge.txt"

# create directories to store outputs from plink and gcta
mkdir Plink_files
plink_out="lambs"

mkdir QC
qc_out="lambs_afterQC"


#Then run plink to merge all files together;
plink --merge-list $output_file  --make-bed --allow-extra-chr --sheep --out Plink_files/$plink_out --allow-no-sex
# founndout that .fam has some IIds like UK17906912232. adding 0 in between
awk '{
    # Detect if IID starts with UK1790691 but missing the 0 after it
    if ($2 ~ /^UK1790691[1-9]/) {
        $2 = substr($2, 1, 9) "0" substr($2, 10)
    }
    print
}' Plink_files/$plink_out.fam > Plink_files/$plink_out"_fixed".fam

# Gnenerate new bed bim files
plink --bfile Plink_files/$plink_out --fam  Plink_files/$plink_out"_fixed".fam --make-bed --out Plink_files/$plink_out"_fixed" --sheep --allow-extra-chr

#create .ped and .map files for the .bed, .bim and .fam files. You need to provide the phenotype file here.
plink --bfile Plink_files/$plink_out"_fixed" --sheep --allow-no-sex --allow-extra-chr --recode --pheno $phenotype --out Plink_files/$plink_out"_fixed"

# get the total number of SNPs in each of the .bim file.
for i in Plink_files/$plink_out"_fixed"*.bim; do cat ${i} | wc -l; done

# identify chromosomes in bim file. below code will list all chromosome codes in your data.
awk '{print $1}' Plink_files/$plink_out"_fixed"*.bim | sort | uniq -c

# create a file with autosomes only (1 to 26)
seq 1 26 > autosomes.txt

# extract autosomal variants only
plink  --bfile Plink_files/$plink_out"_fixed" \
	--chr-set 26 \
	--autosome \
	--make-bed \
	--out Plink_files/$plink_out"_autosomes_only" \
	--sheep \
	--allow-extra-chr

# some preliminary checks on plink files
plink --bfile Plink_files/$plink_out"_fixed"  --not-chr 0,OAR,CONTIG --make-bed --out Plink_files/$plink_out"_CHRfiltered"  --allow-extra-chr --sheep --allow-no-sex
plink --bfile Plink_files/$plink_out"_CHRfiltered"	--not-chr CONTIG --make-bed --out Plink_files/$plink_out"_CHRfiltered"  --allow-extra-chr --sheep --allow-no-sex

plink --bfile Plink_files/$plink_out"_CHRfiltered"  --missing --out Plink_files/$plink_out"_miss_report" --allow-extra-chr --sheep --allow-no-sex
plink --bfile Plink_files/$plink_out"_CHRfiltered"  --het --out Plink_files/$plink_out"_het_report" --allow-extra-chr --sheep --allow-no-sex
plink --bfile Plink_files/$plink_out"_CHRfiltered"  --check-sex --out Plink_files/$plink_out"_sex_report" --allow-extra-chr --sheep --allow-no-sex

# generate QC and summary stats
awk '{print $1, $2, $6}' Plink_files/Birth_weight_miss_report.imiss > Plink_files/fail_reasons.txt
awk '{print $1, $2, $3}' Plink_files/Birth_weight_het_report.het >> Plink_files/fail_reasons.txt
grep "PROBLEM" Plink_files/Birth_weight_sex_report.sexcheck >> Plink_files/fail_reasons.txt


#Next, we run the PLINK command in the terminal. The command that was run is shown below:
#1. Basically, we merge multiple .bed .bim and .fam files for all the genotypes into a single file
#2. Then run a QC on the merged genotype file
#3. get the number of SNPs removed in each filtering step, such as number of snps removed by --maf 0.05 etc..

plink --bfile Plink_files/$plink_out"_CHRfiltered" \
	--maf 0.05 \
	--geno 0.1 \
	--mind 0.1 \
	--hwe 1e-06 \
	--indep-pairwise 50 5 0.8 \
	--allow-extra-chr \
	--allow-no-sex \
	--make-bed \
	--sheep \
	--out QC/$qc_out 
	
# If running pruning in the previous step, then extract the pruned SNPs
plink --bfile QC/$qc_out  --extract QC/$qc_out".prune.in" --out QC/$qc_out"_pruned" --make-bed --sheep --allow-extra-chr

#4. Additional QC that could be run are heterozygosity or inbreeding coefficient. However, this is not necessary for sheep/animal samples as inbreeding is normal.
plink --bfile QC/$qc_out"_pruned"  --extract QC/$qc_out".prune.in" --het  --out QC/$qc_out --sheep --allow-extra-chr

awk -F ' ' '{print $2}' QC/lambs_afterQC_pruned.fam > QC/IDs_passed_PLINK_QC.txt
