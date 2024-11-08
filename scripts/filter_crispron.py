import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("--path")
parser.add_argument("--distance", type=int)
parser.add_argument("--score_threshold", type=float)

args = parser.parse_args()

guides_file = f"{args.path}/output_guides.{args.distance}bp.bed"
output_file = f"{args.path}/output_guides.{args.distance}bp.filtered.bed"

crispron_table = pd.read_table(guides_file, sep="\t", names=['chr','start','end','sequence','score','strand'])

crispron_table = crispron_table[crispron_table.score >= args.score_threshold]

crispron_table.to_csv(output_file, sep="\t", index=False, header=False)