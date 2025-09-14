import subprocess
import os
import pandas as pd
import argparse
import glob


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--log_dir", required=True, help="Directory containing *.log.txt files for all subjects.")
    parser.add_argument("-o", required=True, help="Output file.")
    args = parser.parse_args()

    dict_stats = {}
    indiv = []

    log_files = sorted(glob.glob(os.path.join(args.log_dir, '*.log.txt')))
    for i, filename in enumerate(log_files, 1):
        print(f"\rParsing log {i}/{len(log_files)}", end='' if i < len(log_files) else None)
        participant_id = os.path.split(filename)[1].split('.')[0]
        indiv.append(participant_id)
        if participant_id not in dict_stats:
            dict_stats[participant_id] = {"phased":0, "total":0, "corrected":0}
        with open(os.path.join(args.log_dir, filename), "r") as stream_in:
            for line in stream_in:
                line = line.rstrip()
                #     PHASED  14506 of 151152 all variants (= 0.095970) with at least one other variant
                #     GENOME WIDE PHASE CORRECTED  15 of 151152 variants (= 0.000099)
                xcols = line.split(" ")
                if line.startswith("     PHASED"):
                    dict_stats[participant_id]['phased'] += int(xcols[7])
                    dict_stats[participant_id]['total'] += int(xcols[9])
                elif line.startswith("     GENOME WIDE PHASE CORRECTED"):
                    dict_stats[participant_id]['corrected'] += int(xcols[10])

    df_out = pd.DataFrame(dict_stats.values())

    df_out['p_phased'] = df_out['phased'] / df_out['total']
    df_out['p_corrected'] = df_out['corrected'] / df_out['phased']
    df_out.index = indiv
    df_out.index.name = 'participant_id'
    df_out.to_csv(args.o, sep="\t", index=True, float_format='%.6g')


if __name__ == "__main__":
    main()
