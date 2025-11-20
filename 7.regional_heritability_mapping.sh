# This code assumes that previous plink_gcta script is run.
# Some of the variables from the plink_gcta.sh is used in this script. So make sure that these variables are chosen correctly.

# Include COVARIATES as needed
# If covariates are there, at the end of gcta add the parameter --covar $COVAR_FILE

#!/bin/bash


WORKDIR="/hdd/8TBSSD/Steffi/Sheep_genotypes/all_IDs_2583/ssGWAS/Age_at_slaughter"
cd $WORKDIR/ || exit

OUTDIR="RHM_results"
mkdir $OUTDIR

rhm_phenotype="phenotypes.txt"
COVAR_FILE="covariates.txt"

qc_out="QC/228_UKIDs_afterQC_pruned"

for phenotype in "Age_at_slaughter" 
do
echo "mkdir -p ./RHM_results/${phenotype}"
mkdir -p ./RHM_results/${phenotype}

for win in 500 100 50 20 10
do
step=$((win/2))
echo $win $step
for Ch in {1..26}
do
nsnp=$(wc -l SNP.ch${Ch} | awk '{total += $1} END{print total}')
nwin=$((nsnp/step))
c=0
echo "nsnps" $nsnp
echo "nwin" $nwin
while [ $c -le $nwin ]
do
c=`expr $c + 1`
start=$(($c -1))
start=$(($start * $step))
start=$(($start + 1))
end=$(($start + $win))
end=$(($end - 1))
if [ -f ./RHM_results/${phenotype}/${phenotype}${Ch}_${start}_${end}.hsq ]
then  
echo "nwin" ${phenotype} $c "start" $start "end" $end "exist"
else
echo "win" ${phenotype} $c "start" $start "end" $end
head -n${end} SNP.ch${Ch} | tail -n${win} >> SNP.list${Ch}_${start}_${end}
gcta64 \
--bfile  $qc_out"_phenoAdded" \
--make-grm \
--out ${phenotype}${Ch}_${start}_${end} \
--extract SNP.list${Ch}_${start}_${end} \
--autosome-num 26
echo "Rgrm was calculated"
echo ${phenotype}${Ch}_${start}_${end} > ${phenotype}${Ch}_${start}_${end}.txt
echo "$qc_out""_GRM_autosomes" >>${phenotype}${Ch}_${start}_${end}.txt
echo "${phenotype}${Ch}_${start}_${end}.txt"
head ${phenotype}${Ch}_${start}_${end}.txt
####analysis TB
gcta64 \
--reml \
--mgrm ${phenotype}${Ch}_${start}_${end}.txt \
--reml-lrt 1 \
--pheno $rhm_phenotype \
--covar "$COVAR_FILE" \
--out ./RHM_results/${phenotype}/${phenotype}${Ch}_${start}_${end} 
rm ${phenotype}${Ch}_${start}_${end}.*

fi
done
done
done
done



# Exrtract RHM results
WORKDIR="/hdd/8TBSSD/Steffi/Sheep_genotypes/all_IDs_2583/ssGWAS/Age_at_slaughter"
OUTDIR="RHM_results"
cp SNP* $OUTDIR
cd "$WORKDIR"

for phenotype in "Age_at_slaughter"
do
  echo $phenotype

  # Define window sizes
  for winsize in 500 100 50 20 10
  do
    winover=$((winsize / 2))
    result_file="${phenotype}_RHM_${winsize}_${winover}.Results"

    # Write header
    echo -e "Chr\tStart\tEnd\tV_G1\tV_G2\tV_e\tVp\tH2\tLogL\tlogL0\tLRT\tPval\tN" > "$result_file"
    echo "Processing window size: $winsize, overlap: $winover"

    for Ch in {1..26}
    do
      start=1
      NSNPS=$(wc -l < ./SNP.ch${Ch})
      wsize=$((winsize - winover))
      end=$((NSNPS / wsize))

      for i in $(seq $start $end)
      do
        s0=$((i + 1))
        s1=$((s0 * winover))
        swin=$((s1 + 1 - winsize))
        ewin=$((swin + winsize - 1))

        echo -n "$Ch$(printf '\t')$swin$(printf '\t')$ewin" >> "$result_file"

        hsq_file="./RHM_results/${phenotype}/${phenotype}${Ch}_${swin}_${ewin}.hsq"

        if [ -f "$hsq_file" ]; then
          echo "Found: $hsq_file"

          cnt=0
          while IFS= read -r line
          do
            cnt=$((cnt + 1))
            a=($line)

            case $cnt in
              2) echo -n "$(printf '\t')${a[1]}" >> "$result_file" ;; # V(G1)
              3) echo -n "$(printf '\t')${a[1]}" >> "$result_file" ;; # V(G2)
              4) echo -n "$(printf '\t')${a[1]}" >> "$result_file" ;; # V(e)
              5) echo -n "$(printf '\t')${a[1]}" >> "$result_file" ;; # Vp
              9) echo -n "$(printf '\t')${a[3]}" >> "$result_file" ;; # V(G)/Vp (sum)
              10) echo -n "$(printf '\t')${a[1]}" >> "$result_file" ;; # logL
              11) echo -n "$(printf '\t')${a[1]}" >> "$result_file" ;; # logL0
              12) echo -n "$(printf '\t')${a[1]}" >> "$result_file" ;; # LRT
              14) echo -n "$(printf '\t')${a[1]}" >> "$result_file" ;; # Pval
              15) echo -n "$(printf '\t')${a[1]}" >> "$result_file" ;; # n
            esac
          done < "$hsq_file"

          echo "" >> "$result_file"
        else
          echo "Missing: $hsq_file"
          echo "" >> "$result_file"
        fi
      done
    done
  done
done

# Results will be in workdir/
