# Run these codes in R prior to running this scriptR

WORKDIR="/hdd/8TBSSD/Steffi/Sheep_genotypes/all_IDs_2583/ssGWAS/Age_at_slaughter"
cd $WORKDIR

read.table("QC/lambs_afterQC_phenoAdded.map")->map
for ( i in 1:26){write.table(map[which(map$V1==i),2], file=paste0("SNP.ch",i),sep=" ",quot=F,col.names=F,row.names=F)}
