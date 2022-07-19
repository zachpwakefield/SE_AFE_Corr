README

Within:

1. BAM download using sratoolkit, rMATS calculation, HIT Index calculation

[Outputs RNA-Seq BAM files, rMATS output, HIT Index output.]

a. Change -P, -N flags to project name and BU SCC run name
b. Change initial cd call within run_pipe.sh to the directory that finPipe is in.
c. Change path to --gap_token to the path for your token in gap_fetch_bam call
d. Edit selectSRA.txt to include a list (separated by newlines) of the SRR codes for the RNA-Seq samples
e. Alter any parameters in HIT Index to change in HIT__identity_parameters.txt 
d. To run pipeline (from finPipe/): qsub run_pipe.sh

(Option to run with FASTQ or run Kallisto not developed)


2. SE & AFE Kendall tau correlation calculations

[Results in various outputs: seBED and afeBED files for each sextile. Distribution across KB correlation value, data frame with tau for each se & are pair]

a. Edit /finPipe/Corr/parPipe.sh with correct sample and output directory (output should be /finPipe/Corr/corrOut/ and sample directory /finPipe/runs/)
b. To run pipeline (from finPipe/Corr/): 