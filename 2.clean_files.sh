#In the terminal:
WORKDIR="/hdd/8TBSSD/Steffi/Sheep_genotypes/all_IDs_2583/ssGWAS/Age_at_slaughter"
cd $WORKDIR/

#if you have multiple zip files
for i in *.zip; do unzip ${i}; done

#Then remove .log and .nosex files
rm *.zip *.log *.nosex

#if you have .txt in file names; remove them
rename 's/\.txt//' *.txt.*
