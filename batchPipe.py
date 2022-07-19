import argparse
import glob
import os.path
import re
import subprocess


def gap_fetch_fastq(runs, gap_token):
    with open(runs, 'r') as fh:
        runs = [line.strip('\n') for line in fh.readlines()]

    for run in runs:
        cmd = (
            f"prefetch --ngc {gap_token} {run} --max-size 30g; "
            f"cd {run}; "
            f"fasterq-dump --ngc {gap_token} {run}.sra; "
            f"cd ..")

        subprocess.call(cmd, shell=True)

        # Remove extraneous data.
        for path in glob.iglob(f"{run}/*"):
            ext = re.search(r'.\w*$', path).group()
            if ext != ".fastq":
                cmd = f"rm -rf {path}"
                subprocess.call(cmd, shell=True)


def gap_fetch_bam(runs, gap_token):
    with open(runs, 'r') as fh:
        runs = [line.strip('\n') for line in fh.readlines()]

    for run in runs:
        cmd = (
            f"prefetch --ngc {gap_token} {run} --max-size 40g; "
            f"cd {run}; "
            f"sam-dump --ngc {gap_token} {run}.sra | "
            f"samtools view -b - > {run}.bam; "
            f"cd ..")

        subprocess.call(cmd, shell=True)

        # Remove extraneous data.
        for path in glob.iglob(f"{run}/*"):
            ext = re.search(r'.\w*$', path).group()
            if ext != ".bam":
                cmd = f"rm -rf {path}"
                subprocess.call(cmd, shell=True)


def sort_index(run_dir):
    for path in glob.iglob(f'{run_dir}/*'):
        run = os.path.basename(path)
        cmd = (
            f"samtools sort {path}/{run}.bam -o {path}/{run}_sorted.bam; "
            f"samtools index {path}/{run}_sorted.bam; ")

        subprocess.call(cmd, shell=True)


def kallisto(run_dir, idx):
    for path in glob.iglob(f"{run_dir}/*"):
        subprocess.call(
            f"kallisto quant -i {idx} -o {path}/kallisto {path}/fastq/*.fastq",
            shell=True)


def rmats(run_dir, gtf):
    for path in glob.iglob(f"{run_dir}/*"):
        run = os.path.basename(path)
        cmd = (
            f"mkdir {path}/rmats {path}/rmats/tmp; "
            f"echo '{path}/{run}.bam' > {path}/bam_path.txt; "
            f"rmats.py --b1 {path}/bam_path.txt --gtf {gtf} -t paired "
            f"--readLength 76 --nthread 4 --od {path}/rmats "
            f"--tmp {path}/rmats/tmp --statoff; "
            f"rm -rf {path}/bam_path.txt {path}/rmats/tmp")

        subprocess.call(cmd, shell=True)


def hit_index(run_dir, bed, bootstrap):
    for path in glob.iglob(f"{run_dir}/*"):
        run = os.path.basename(path)

        cmd = (
            f"python HITindex_classify.py "
            f"--junctionReads --bam {path}/{run}.bam "
            f"--juncbam {path}/{run}.juncs.bam "
            f"--readstrand fr-unstrand "
            f"--bed {bed} "
            f"--outname {path}/{run} "
            f"--overlap 15 "
            f"--readnum 10 "
            f"--bootstrap {bootstrap} "
            f"--classify "
            f"--HITindex")

        subprocess.call(cmd, shell=True)


def calculate_PSI(run_dir, bed):
    for path in glob.iglob(f"{run_dir}/*"):
        run = os.path.basename(path)

        cmd = (
            f"python HITindex_classify.py "
            f"--calculatePSI "
            f"--metricsID {path}/hit/{run}.exon "
            f"--outname {path}/hit/{run}")

        subprocess.call(cmd, shell=True)




if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    # SAM/FASTQ fetch arguments.
    parser.add_argument('--function')
    parser.add_argument('--gap_token')
    parser.add_argument('--runs')

    # Tool arguments.
    parser.add_argument('--run_dir')
    parser.add_argument('--idx')
    parser.add_argument('--gtf')
    parser.add_argument('--bed')
    parser.add_argument('--bootstrap')

    args = parser.parse_args()

    if args.function == 'gap_fetch_fastq':
        gap_fetch_fastq(args.runs, args.gap_token)

    elif args.function == 'gap_fetch_bam':
        gap_fetch_bam(args.runs, args.gap_token)

    elif args.function == 'sort_index':
        sort_index(args.run_dir)

    elif args.function == 'kallisto':
        kallisto(args.run_dir, args.idx)

    elif args.function == 'rmats':
        rmats(args.run_dir, args.gtf)

    elif args.function == 'hit_index':
        hit_index(args.run_dir, args.bed, args.bootstrap)

    elif args.function == 'calculate_PSI':
        calculate_PSI(args.run_dir, args.bed)

    else:
        print(
            "--function must be passed either `gap_fetch_fastq`, "
            "`gap_fetch_bam`, `sort_index`, `kallisto`, `hit_index`, `calculate_PSI`, `rmats`.")
