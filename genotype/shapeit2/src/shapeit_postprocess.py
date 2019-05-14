import numpy as np
import argparse
import subprocess
import gzip
from datetime import datetime


def gt_dosage(gt):
    """Convert unphased genotype to dosage"""
    x = gt.split(b'/')
    return int(x[0])+int(x[1])


def pgt_dosage(gt):
    """Convert phased genotype to dosage"""
    x = gt.split(b'|')
    return int(x[0])+int(x[1])


def calculate_missingness(miss_buf):
    if len(miss_buf)==1:
        pct_missing = sum(miss_buf[0]) / len(miss_buf[0])
    else:
        pct_missing = np.all(miss_buf, axis=0).sum() / len(miss_buf[0])
    return pct_missing


def patch_phased_vcf(vcf_file, phased_vcf_file, patched_vcf_file, missingness_threshold=0.15):
    """
    Patches phased VCF produced by SHAPEIT2
    - if the dosage is different between the phased and unphased VCF, assigns to missing
    - for split biallelic sites, assigns to missing if ALT site is different (e.g., 0|2)
    - filters out sites with missingness > missingness_threshold
    """
    assert patched_vcf_file.endswith('.vcf.gz')
    log_file = patched_vcf_file.replace('.vcf.gz', '.log')

    pgt_set = set([b'.|.', b'0|0', b'0|1', b'1|0', b'1|1'])

    num_var = int(subprocess.check_output('zcat {} | grep -v "#" | wc -l'.format(phased_vcf_file), shell=True).decode())
    print('  * parsing {} sites'.format(num_var))

    bgzip = subprocess.Popen('bgzip -c > '+patched_vcf_file, stdin=subprocess.PIPE, shell=True)
    nwritten = 0
    ndropped = 0
    with gzip.open(vcf_file) as vcf, gzip.open(phased_vcf_file) as phased_vcf, open(log_file, 'w') as log:
        # skip header of unphased VCF
        for line in vcf:
            if line[:6]==b'#CHROM':
                break
        # copy header of phased VCF
        for pline in phased_vcf:
            if pline[:6]==b'#CHROM':
                bgzip.stdin.write(b'##Note=Processed using shapeit_postprocess.py\n')
                bgzip.stdin.write(pline)
                break
            bgzip.stdin.write(pline)
        assert line==pline

        # iterate through variants
        buf = []
        miss_buf = []
        previous_pos = None
        for k, (line, pline) in enumerate(zip(vcf, phased_vcf), 1):
            line = line.strip().split(b'\t')
            s = pline.strip().split(b'\t')
            assert line[2]==s[2]  # same variant

            pos = line[0]+b'_'+line[1]

            # if new site, write previous
            if pos != previous_pos and previous_pos is not None:
                # calculate missingness and write site
                pct_missing = calculate_missingness(miss_buf)
                if pct_missing <= missingness_threshold:
                    for b in buf:
                        bgzip.stdin.write(b'\t'.join(b)+b'\n')
                    nwritten += len(buf)
                else:
                    log.write('{}: missingness {:.4f}, removed (n = {}).\n'.format(previous_pos.decode(), pct_missing, len(buf)))
                    ndropped += len(buf)

                # reset buffers
                buf = []
                miss_buf = []

            # parse and process current site
            gt = line[9:]
            pgt = s[9:]

            # check for dosage differences at non-imputed sites, set to missing
            ix = [i!=b'./.' and gt_dosage(i)!=pgt_dosage(j) for i,j in zip(gt, pgt)]
            if any(ix):  # dosage mismatch, set to missing
                pgt = [b'.|.' if i else g for i,g in zip(ix,pgt)]
                log.write('{}: mismatched dosage for {} samples (set to missing).\n'.format(s[2].decode(), sum(ix)))

            # reset split sites
            if b'wasSplit' in line[7]:
                # split multi-allelic sites: a sample can only have genotype at one split site; others must be set to missing
                # Since all are imputed, set back to missing (if all were originally missing, imputed call is ambiguous).
                ix = [g==b'./.' for g in gt]
                if any(ix):
                    pgt = [b'.|.' if i else g for i,g in zip(ix,pgt)]
                    log.write('{}: split site; {} imputed sites reverted to missing.\n'.format(s[2].decode(), sum(ix)))

            buf.append(s[:9]+pgt)
            miss_buf.append([i==b'./.' or j==b'.|.' for i,j in zip(gt, pgt)])
            previous_pos = pos

            if np.mod(k, 10000)==0:
                print('\r    * variants processed: {}'.format(k), end='', flush=True)
        print('\r    * variants processed: {}'.format(k), end='', flush=True)
        print()

        # last site: calculate missingness and write site
        pct_missing = calculate_missingness(miss_buf)
        if pct_missing <= missingness_threshold:
            for b in buf:
                bgzip.stdin.write(b'\t'.join(b)+b'\n')
            nwritten += len(buf)
        else:
            log.write('{}: missingness {:.4f}, removed (n = {}).\n'.format(previous_pos.decode(), pct_missing, len(buf)))
            ndropped += len(buf)

    stdout, stderr = bgzip.communicate()
    print('  * wrote {} sites'.format(nwritten))
    print('  * dropped {} sites'.format(ndropped))


if __name__=='__main__':

    parser = argparse.ArgumentParser(description='Run post-processing for SHAPEIT2 phasing')
    parser.add_argument('vcf_file', type=str, help='Unphased VCF')
    parser.add_argument('phased_vcf_file', type=str, help='SHAPEIT2 output')
    parser.add_argument('patched_vcf_file', help='Output VCF')
    parser.add_argument('--missingness_threshold', default=0.15, type=np.float64, help='Missingness threshold')
    args = parser.parse_args()

    print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Patching VCF', flush=True)
    patch_phased_vcf(args.vcf_file, args.phased_vcf_file, args.patched_vcf_file, missingness_threshold=args.missingness_threshold)

    # index
    print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Indexing patched VCF', flush=True)
    subprocess.check_call('tabix '+args.patched_vcf_file, shell=True)
    n = int(subprocess.check_output('bcftools index -n {}'.format(args.patched_vcf_file), shell=True).decode().strip())
    print('  * {} sites remaining'.format(n))
