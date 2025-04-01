#!/bin/bash

# %%%%%%%%%%%%%%%%%%%%%%%%
# %%% PICO: RUN PLINK2 %%%
# %%%%%%%%%%%%%%%%%%%%%%%%

# This script is used for the Vitis labrusca dataset (929 samples, 32,499 SNPs) obtained from Dr. Zoe Migicovsky and the USDA. It coverts the .ped and .map file
# to a VCF file (using the --recode-structure argument), which is then converted to a genind file in R.

# Arguments:
# --pedmap: specifies the prefixes of the .ped and .map files to pass to PLINK2. This files were obtained directly from Dr. Migicovsky.
#
# --recode: specifies the format of the output (VCF)
#
# --out: specify name of output VCF file

plink2 --pedmap pcjl_grape_maf001_min2alleles_noindels_filtered_impute_hwe_lab --recode vcf --out VILA

