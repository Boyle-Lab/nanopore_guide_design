import argparse
import glob
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("--path")
parser.add_argument("--distance")
parser.add_argument("--min_gc", type=int)
parser.add_argument("--max_gc", type=int)
parser.add_argument("--max_self_complementarity", type=int)
parser.add_argument("--min_efficiency_score", type=float)
parser.add_argument("--max_mm0", type=int)
parser.add_argument("--max_mm1", type=int)
parser.add_argument("--max_mm2", type=int)
parser.add_argument("--max_mm3", type=int)

args = parser.parse_args()

distance = args.distance
min_gc=args.min_gc
max_gc=args.max_gc
max_self_complementarity=args.max_self_complementarity
min_efficiency_score=args.min_efficiency_score
max_mm0=args.max_mm0
max_mm1=args.max_mm1
max_mm2=args.max_mm2
max_mm3=args.max_mm3

output_path=args.path
output_files=glob.glob(f'{output_path}/*.txt')

for output_file in output_files:
    name = output_file.split('.txt')[0]
    filtered_file = f'{name}.filtered.tsv'
    bedfile = f'{name}.bed'
    output_table = pd.read_table(output_file, sep="\t", names=['ID','Sequence','Location','Strand','GC','SelfComp','mm0','mm1','mm2','mm3','EfficiencyScore'])
    ## Filter CHOPCHOP output
    #Filter on GC
    output_table = output_table[(output_table['GC'] >= min_gc) & (output_table['GC'] <= max_gc)]
    #
    #Filter on SelfComp
    output_table = output_table[(output_table['SelfComp']) <= max_self_complementarity]
    #
    #Filter on mm0/mm1/mm2/mm3
    try:
        output_table['mm0'] = output_table['mm0'].str.replace(">=","")
    except:
        pass

    try:
        output_table['mm1'] = output_table['mm1'].str.replace(">=","")
    except:
        pass

    try:
        output_table['mm2'] = output_table['mm2'].str.replace(">=","")
    except:
        pass

    try:
        output_table['mm3'] = output_table['mm3'].str.replace(">=","")
    except:
        pass

    #
    output_table['mm0'] = output_table['mm0'].astype('int64')
    output_table['mm1'] = output_table['mm1'].astype('int64')
    output_table['mm2'] = output_table['mm2'].astype('int64')
    output_table['mm3'] = output_table['mm3'].astype('int64')
    #
    output_table = output_table[(output_table['mm0'] <= max_mm0) & (output_table['mm1'] <= max_mm1) &
                                (output_table['mm2'] <= max_mm2) & (output_table['mm3'] <= max_mm3)]
    #
    #Write filtered table
    if output_table.size == 0:
        continue
    output_table.to_csv(filtered_file, sep="\t", index=False)
    #
    # Build BED files
    output_table[['chr','position']] = output_table['Location'].str.split(':', n=1, expand=True)
    output_table['position'] = output_table['position'].astype('int64')
    output_table['name'] = '.'
    output_table['score'] = '.'
    #
    plus_strand = (output_table.Strand == "+")
    minus_strand = (output_table.Strand == "-")
    #
    output_table.loc[plus_strand,'mer_start'] = output_table['position'] - 5
    output_table.loc[plus_strand,'mer_end'] = output_table['position'] + 25
    output_table.loc[minus_strand,'mer_start'] = output_table['position'] - 4
    output_table.loc[minus_strand,'mer_end'] = output_table['position'] + 26
    #
    output_table['mer_start'] = output_table['mer_start'].astype('int64')
    output_table['mer_end'] = output_table['mer_end'].astype('int64')
    #
    bed_table = output_table[['chr','mer_start','mer_end','name','score','Strand']].copy()
    #
    bed_table.to_csv(bedfile, sep="\t", index=False, header=False)