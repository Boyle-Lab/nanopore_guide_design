import argparse
import glob
import numpy as np
import os
import pandas as pd

def read_targets_bed(targets_file):
    targets = pd.read_table(targets_file, header=None,names=['chr','start','end'])
    targets['region'] = targets['chr'] + ":" + targets['start'].astype('str') + "-" + targets['end'].astype('str')
    targets.set_index('region',inplace=True)
    return(targets)

def build_targets_guides(targets):
    targets_guides = targets.copy()

    targets_guides['U_chr'] = ""
    targets_guides['U_start'] = np.nan
    targets_guides['U_end'] = np.nan
    targets_guides['U_seq'] = ""
    targets_guides['U_score'] = np.nan
    targets_guides['U_strand'] = ""
    targets_guides['U_distance'] = np.nan

    targets_guides['D_chr'] = ""
    targets_guides['D_start'] = np.nan
    targets_guides['D_end'] = np.nan
    targets_guides['D_seq'] = ""
    targets_guides['D_score'] = np.nan
    targets_guides['D_strand'] = ""
    targets_guides['D_distance'] = np.nan

    return(targets_guides)

def write_targets_guides(targets_guides, targets_guides_file):
    targets_guides.to_csv(targets_guides_file, sep="\t", index=True, header=True)

# Parameters
parser = argparse.ArgumentParser()
parser.add_argument("--targets-bed")
parser.add_argument("--distance", type=int)
parser.add_argument("--score-threshold", type=float)
parser.add_argument("--output")

args = parser.parse_args()

targets_file = args.targets_bed
distance = args.distance
score_threshold = args.score_threshold
output_path = args.output


# Read in targets BED
targets = read_targets_bed(targets_file)

# Read/build targets and candidate guides file
targets_guides_file = f'{output_path}/targets_candidate_guides.tsv'

if not os.path.isfile(targets_guides_file):
    targets_guides = build_targets_guides(targets)

else:
    targets_guides = pd.read_table(targets_guides_file)
    targets_guides.set_index('region',inplace=True)

# Parse upstream candidate guides
for index, row in targets.iterrows():
    try:
        region = f"{row['chr']}:{row['start']}-{row['end']}"
        upstream_guide_filepath = f"{row['chr']}.{row['start']}.{row['end']}.upstream_high_scoring_guides.{distance}bp.bed"
        upstream_guide_file = pd.read_table(f"{output_path}/{upstream_guide_filepath}", header=None, names=['chr','start','end','seq','score','strand'])
    except:
        continue
    for index, row in upstream_guide_file.iterrows():
        if row['score'] >= score_threshold:
            if targets_guides.loc[[region]]['U_score'].isna().values[0]:
                targets_guides.at[region, 'U_chr'] = row['chr']
                targets_guides.at[region, 'U_start'] = row['start']
                targets_guides.at[region, 'U_end'] = row['end']
                targets_guides.at[region, 'U_seq'] = row['seq']
                targets_guides.at[region, 'U_score'] = row['score']
                targets_guides.at[region, 'U_strand'] = row['strand']
                targets_guides.at[region, 'U_distance'] = distance
            elif targets_guides.loc[[region]]['U_score'].values[0] < row['score']:
                targets_guides.at[region, 'U_chr'] = row['chr']
                targets_guides.at[region, 'U_start'] = row['start']
                targets_guides.at[region, 'U_end'] = row['end']
                targets_guides.at[region, 'U_seq'] = row['seq']
                targets_guides.at[region, 'U_score'] = row['score']
                targets_guides.at[region, 'U_strand'] = row['strand']
                targets_guides.at[region, 'U_distance'] = distance

# Parse downstream candidate guides
for index, row in targets.iterrows():
        try:
            region = f"{row['chr']}:{row['start']}-{row['end']}"
            downstream_guide_filepath = f"{row['chr']}.{row['start']}.{row['end']}.downstream_high_scoring_guides.{distance}bp.bed"
            downstream_guide_file = pd.read_table(f"{output_path}/{downstream_guide_filepath}", header=None, names=['chr','start','end','seq','score','strand'])
        except:
            continue
        for index, row in downstream_guide_file.iterrows():
            if row['score'] >= score_threshold:
                if targets_guides.loc[[region]]['D_score'].isna().values[0]:
                    targets_guides.at[region, 'D_chr'] = row['chr']
                    targets_guides.at[region, 'D_start'] = row['start']
                    targets_guides.at[region, 'D_end'] = row['end']
                    targets_guides.at[region, 'D_seq'] = row['seq']
                    targets_guides.at[region, 'D_score'] = row['score']
                    targets_guides.at[region, 'D_strand'] = row['strand']
                    targets_guides.at[region, 'D_distance'] = distance
                elif targets_guides.loc[[region]]['D_score'].values[0] < row['score']:
                    targets_guides.at[region, 'D_chr'] = row['chr']
                    targets_guides.at[region, 'D_start'] = row['start']
                    targets_guides.at[region, 'D_end'] = row['end']
                    targets_guides.at[region, 'D_seq'] = row['seq']
                    targets_guides.at[region, 'D_score'] = row['score']
                    targets_guides.at[region, 'D_strand'] = row['strand']
                    targets_guides.at[region, 'D_distance'] = distance

# Write output
write_targets_guides(targets_guides, targets_guides_file)