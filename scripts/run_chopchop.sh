#!/bin/bash

DISTANCE="$2"
FLANK="$3"
OUTPUT="$4"
CHOPCHOP="$5"

if [ ! -d "$OUTPUT" ]; then
    mkdir -p "$OUTPUT"
fi

# # loop through each line in the file
IFS=$'\t' read -r col1 col2 col3 <<< "$1"
# make the second and third columns integer variables, if possible
if [[ $col2 =~ ^[0-9]+$ ]]; then
  col2_int=$((col2))
else
  col2_int="not an integer"
  echo 'stop coord "$col2" is "$col2_int"'
  exit
fi
col3="$(echo -e "${col3}" | sed -e 's/[[:space:]]*$//')"
if [[ $col3 =~ ^[0-9]+$ ]]; then
  col3_int=$((col3))
else
  col3_int="not an integer"
  echo 'stop coord "$col3" is "$col3_int'
  exit
fi
############################## Makes regions and pass to chopchop

upstream_coods_start=$(($col2_int - $DISTANCE - $FLANK))
upstream_coods_end=$(($col2_int - $DISTANCE + $FLANK))

downstream_coods_start=$(($col3_int + $DISTANCE-$FLANK))
downstream_coods_end=$(($col3_int + $DISTANCE+$FLANK))

LOC="$col1"."$col2_int"."$col3_int"

# echo "$col1":"$upstream_coods_start"-"$upstream_coods_end"
# echo "$col1":"$downstream_coods_start"-"$downstream_coods_end"

#Upstream
python "$CHOPCHOP" \
  -Target "$col1":"$upstream_coods_start"-"$upstream_coods_end" \
  -BED \
  -GenBank \
  -G hg38 \
  -filterGCmin 20 \
  -filterGCmax 80 \
  -filterSelfCompMax 0 \
  -t WHOLE \
  -n N \
  -R 4 \
  -T 1 \
  -g 20 \
  -scoringMethod DOENCH_2016 \
  -f NN \
  -v 3 \
  -M NGG \
  -BB AGGCTAGTCCGT \
  -o "$OUTPUT"/"$LOC"_temp \
  --nonOverlapping | awk '$4 == "+" {print}' | head -n 25 >> "$OUTPUT"/"$LOC".upstream.txt

#Downstream
python "$CHOPCHOP" \
  -Target "$col1":"$downstream_coods_start"-"$downstream_coods_end" \
  -BED \
  -GenBank \
  -G hg38 \
  -filterGCmin 20 \
  -filterGCmax 80 \
  -filterSelfCompMax 0 \
  -t WHOLE \
  -n N \
  -R 4 \
  -T 1 \
  -g 20 \
  -scoringMethod DOENCH_2016 \
  -f NN \
  -v 3 \
  -M NGG \
  -BB AGGCTAGTCCGT \
  -o "$OUTPUT"/"$LOC"_temp \
  --nonOverlapping | awk '$4 == "-" {print}' | head -n 25 >> "$OUTPUT"/"$LOC".downstream.txt
rm -r "$OUTPUT"/"$LOC"_temp
