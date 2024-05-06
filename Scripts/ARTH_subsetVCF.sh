#!/bin/bash

# %%%%%%%%%%%%%%%%%%%%%%%%
# %%% ARTH: SUBSET VCF %%%
# %%%%%%%%%%%%%%%%%%%%%%%%

# This script was used for the Arabidopsis thaliana dataset (1,135 samples), sourced from 1001 Genomes Consortium, 2016.
# It creates the ARTH_10k.vcf file which is read in during the geographic-genetic resampling steps in the Demo_ARTH.R script.
# This script subsets the original VCF file containing full genome information (here: https://1001genomes.org/data/GMI-MPI/releases/v3.1/1001genomes_snp-short-indel_only_ACGTN.vcf.gz)
# to just 10,000 SNPs. NOTE: this script is only run once, and is not a part of the geographic-genetic correlation workflow; it's kept here
# for documentation purposes.

# Begin by copying the header information of the original VCF file to a new VCF file
head -n 9 ARTH_fullGenome.vcf > ARTH_10k.vcf

# Subsetthe original VCF, selecting 10,000 random loci, and append the output to the previously created VCF
grep -v "^#" ARTH_fullGenome.vcf | uniq | shuf | head -n 10000 | sort -n >> ARTH_10k.vcf

# A description of the steps involved in this command are provided below

# 1. grep pulls out all the lines in the VCF except the commented header lines (those starting with #)
# 2. uniq reduces them to a unique list of lines (likely unnecessary, but included for safety reasons)
# 3. shuf randomly shuffle those lines
# 4. head is used to take the first 10,000 of the randomly shuffled lines
# 5. sort is used to sort them again
