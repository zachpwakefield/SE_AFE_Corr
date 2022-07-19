#!/bin/bash -l


# Add list of SRR's and ensure the gtf files are in correct format for both rMATS and HIT Index

#$ -P evolution
#$ -l h_rt=99:59:00
#$ -N fullPipe
#$ -j y
#$ -m ea
#$ -pe omp 28
#$ -l mem_per_core=18G


# Keep track of information related to the current job
echo "=========================================================="
echo "Start date : $(date)"
echo "Job name : $JOB_NAME"
echo "Job ID : $JOB_ID  $SGE_TASK_ID"
echo "=========================================================="

module load samtools
module load sratoolkit
module load miniconda/4.9.2
module load rmats/4.1.1
module load bedtools/2.30.0
cd /projectnb2/evolution/zwakefield/finPipe
cd runs/

# Download BAM
python3 ../batchPipe.py --function gap_fetch_bam --runs ../selectSRA.txt --gap_token /projectnb2/evolution/zwakefield/finPipe/token.ngc

echo "BAM DOWNLOADS COMPLETE"

cd ../

# Run rMATS
python3 batchPipe.py --function rmats --run_dir runs/ --gtf hg19.ensGene.gtf

echo "RMATS COMPLETE"

module unload rmats/4.1.1
module unload miniconda/4.9.2
module load python3

pip install scipy
pip install numpy
pip install pysam
pip install pybedtools
pip install pandas
pip install pymc3

# Generate Meta Exons`
python3 HITindex_annotate.py --gtf hg19.gene.gtf.txt --ss3buffer 50 --ss5buffer 20 --outfile HITmeta/HITmeta

echo "HIT META COMPLETE"

# Run HIT Index
python3 batchPipe.py --function hit_index --run_dir runs/ --bed HITmeta/HITmeta_ss3-50ss5-20.buffer --bootstrap 5

cd runs/
for f in ./*; do mkdir $f/hit; mv $f/*.exon $f/hit/; mv $f/*.AFEPSI $f/hit/; mv $f/*.ALEPSI $f/hit/; mv $f/*.junc* $f/hit/; done
cd ..

python3 batchPipe.py --function calculate_PSI --run_dir runs/ --bed HITmeta/HITmeta_ss3-50ss5-20.buffer


echo "HIT INDEX COMPLETE"


mv *.o* outMessages/
mv *.po* outMessages/
