import subprocess
import os
import glob
import pandas as pd
import multiprocessing as mp
import subprocess
import argparse


def read_result(xinput):
    index, path = xinput
    print(f"  * {index}: {path}")
    df_ind = pd.read_csv(path, sep="\t")
    df_ind.rename(columns={'bam': 'sample_id'}, inplace=True)
    df_ind['sample_id'] = df_ind['sample_id'].apply(lambda x: x.split('.')[0])
    # add individual
    df_ind['individual_id'] = df_ind['sample_id'].apply(
        lambda x: x.split('-')[0] if x.startswith('DGTEX') else '-'.join(x.split('-')[:2]))

    out = []
    for sample_id in df_ind['sample_id'].unique():
        result_all = []
        result_phased = []
        df_sample = df_ind[(df_ind['sample_id'] == sample_id)]

        for gw_phased, aCount, bCount in zip(df_sample['gw_phased'], df_sample['aCount'], df_sample['bCount']):
            result_all.append(f"{aCount}|{bCount}")

            if int(gw_phased) == 1:
                result_phased.append(f"{aCount}|{bCount}")
            else:
                result_phased.append("0|0")

        with open(f"tmp/all/{sample_id}.txt", "w") as f:
            f.write("\n".join([sample_id] + result_all) + "\n")

        with open(f"tmp/gw_phased/{sample_id}.txt", "w") as f:
            f.write("\n".join([sample_id] + result_phased) + "\n")

        out.append(sample_id)

    return out


def parallelize(function, threads, pool_input):
    if len(pool_input) > 0:
        threads = min([len(pool_input),threads])
        if threads > 1:
            pool = mp.Pool(processes=threads)
            pool_output = pool.map(function, pool_input)
            pool.close() # no more tasks
            pool.join()  # wrap up current tasks
        else:
            pool_output = []
            for input in pool_input:
                pool_output.append(function(input))
    else:
        pool_output = []

    return pool_output


def generate_basic_matrix(path):
    xdf = pd.read_csv(path, sep="\t", index_col=False)
    one_bam = list(set(xdf['bam']))[0]
    xdf = xdf[xdf['bam'] == one_bam]
    return xdf[['contig', 'start', 'stop', 'name']].rename(columns={'contig':'#contig'}).reset_index(drop=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", required=True, help="Directory containing *.gene_ae.txt.gz files for all subjects.")
    parser.add_argument("-t", required=False, default=1, type=int, help="Number of threads to use.")
    parser.add_argument("-o", required=True, help="Output file prefix.")
    args = parser.parse_args()

    # make needed temporary directories
    for d in ["tmp", "tmp/all", "tmp/gw_phased"]:
        if not os.path.exists(d):
            os.mkdir(d)

    print("Loading results ...")
    # contig  start   stop    name    aCount  bCount  totalCount      log2_aFC        n_variants      variants        gw_phased       bam
    # 1       11869   14362   ENSG00000223972 0       0       0       inf     0               1       GTEX-YJ89-0011-R9a-SM-4SOK7
    xinput = [[k, f] for k,f in enumerate(sorted(glob.glob(os.path.join(args.i, "*.gene_ae.txt.gz"))), 1)]
    xresults_raw = parallelize(read_result, args.t, xinput)
    xresults = [j for i in xresults_raw for j in i]  # all sample IDs

    header_df = generate_basic_matrix(xinput[0][1])
    print("Saving sample matrix (all)...")
    df_matrix = pd.concat([header_df] + [pd.read_csv(f"tmp/all/{xsample}.txt", sep="\t") for xsample in xresults], axis=1)
    df_matrix.to_csv(f"{args.o}.txt", sep="\t", index=False)
    subprocess.check_call(f"bgzip -f {args.o}.txt", shell=True)
    subprocess.check_call("rm tmp/all/*.txt", shell=True)

    print("Saving sample matrix (gw_phased)...")
    df_matrix = pd.concat([header_df] + [pd.read_csv(f"tmp/gw_phased/{xsample}.txt", sep="\t") for xsample in xresults], axis=1)
    df_matrix.to_csv(f"{args.o}.gw_phased.txt", sep="\t", index=False)
    subprocess.check_call(f"bgzip -f {args.o}.gw_phased.txt", shell=True)
    subprocess.check_call("rm tmp/gw_phased/*.txt", shell=True)



if __name__ == "__main__":
    main()
