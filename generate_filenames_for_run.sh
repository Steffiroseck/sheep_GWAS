#In the terminal:
WORKDIR="/hdd/8TBSSD/Steffi/Sheep_genotypes/all_IDs/"
cd $WORKDIR/

# Next we need the file names for each sample in a format for PLINK to run. Lets do that.
touch $WORKDIR/Filenames_for_plink_merge.txt
output_file="$WORKDIR/Filenames_for_plink_merge.txt"

#Create or clear the output file
> "$output_file"

#Use awk to process file names

ls "$WORKDIR/" | awk -F. '
{
    if ($NF ~ /^(bed|bim|fam)$/) {
        files[$1] = files[$1] ? files[$1] " " $0 : $0
    }
}
END {
    for (id in files) {
        print files[id]
    }
}' | sort > "$output_file"

echo "Grouped filenames saved to $output_file"
