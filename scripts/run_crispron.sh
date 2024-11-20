#!/bin/bash

candidate_bed=$1
OUTPUT=$2
DISTANCE=$3
BEDTOOLS=$4
CRISPRON=$5
name=$(basename $candidate_bed .bed)

# Use bed coords to make fastas for scoring
${BEDTOOLS} getfasta -s -fi /data/genomes/hg38/seq/hg38.fa -bed $candidate_bed > "$OUTPUT"/"$name"."$DISTANCE"bp.fasta

#Score with crisprON TL
${CRISPRON} "$OUTPUT"/"$name"."$DISTANCE"bp.fasta "$OUTPUT"/"$name"

mv "$OUTPUT"/"$name"/crispron.csv "$OUTPUT"/"$name".crispron_results."$DISTANCE"bp.csv
rm -r "$OUTPUT"/"$name"
# Gets the 3 highest scoring guides - no filter for score > 10 
sort -k3 -nr -t, "$OUTPUT"/"$name".crispron_results."$DISTANCE"bp.csv > "$OUTPUT"/"$name".crispron_results."$DISTANCE"bp.sorted.csv
head -n -1 "$OUTPUT"/"$name".crispron_results."$DISTANCE"bp.sorted.csv | awk -F'[(),]|[_]' 'NR < 4 {split($1,a,":|-"); print a[1]"\t"a[2]"\t"a[3]"\t"$6"\t"$7"\t"$2}' > "$OUTPUT"/"$name"_high_scoring_guides."$DISTANCE"bp.bed
