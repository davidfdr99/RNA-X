#!/bin/bash

# Shell script to fetch NGS data given a SRA Accession list with one accession per line

echo "Started downloading SRAs..."
while IFS= read -r 
do 
  prefetch "$line"
done < $1

echo "Finished downloading SRAs, fetching..."

cat $1 | xargs fasterq-dump

while IFS= read -r line
do
    rm -r "$line" &
done < $1

echo "Cleaned."

for f in ./*_{1,2}.fastq; do
    gzip "$f" &
done
