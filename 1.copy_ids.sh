#!/bin/bash

SRC_DIR="/hdd/8TBSSD/Steffi/Sheep_genotypes/afbi_all"
DST_DIR="/hdd/8TBSSD/Steffi/Sheep_genotypes/all_IDs_2583/ssGWAS/Age_at_slaughter/"
FILE_LIST="/hdd/8TBSSD/Steffi/Sheep_genotypes/lambids"

# Create destination directory if it doesn't exist
mkdir -p "$DST_DIR"

# Loop through each filename in the list
while IFS= read -r filename; do
  src_file="$SRC_DIR/$filename"
  if [[ -f "$src_file" ]]; then
    cp "$src_file" "$DST_DIR/"
    echo "Copied: $filename"
  else
    echo "Not found: $filename"
  fi
done < "$FILE_LIST"
