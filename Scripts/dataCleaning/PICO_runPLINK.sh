#!/bin/bash

# %%%%%%%%%%%%%%%%%%%%%%%
# %%% PICO: RUN PLINK %%%
# %%%%%%%%%%%%%%%%%%%%%%%

# This script is used for the Pinus contorta dataset (929 samples, 32,499 SNPs) generated from MacLachlan et al. 2021. It coverts the .ped file
# generated from the supplemental code of the publication to a STRUCTURE file (using the --recode-structure argument), which is then converted to a genind file in R.

# Arguments:
# --file: specifies .ped file to pass to PLINK. This was generated from the first R script shared from MacLachlan et al. 2021, entitled
# 1_MacLachlan_etal_Pine_GPA_ped&mapfile_formatting_Jan10th2021.R
#
# --recode-structure: tells PLINK to rewrite the .ped file as a STRUCTURE input file
#
# --noweb: tells PLINK to skip the web-check (not really sure what the web-check does, but it wasn't able to connect to the Internet on initial runs)
#
# --compound-genotypes: let PLINK know that there's no spacing between alleles
#
# --allow-no-sex: sex of individuals in .ped file isn't specified
#
# --out: specify name of output file

plink --file Pine_NaturalComp_85SNPfilter --recode-structure --noweb --compound-genotypes --allow-no-sex --out Pine_NaturalComp_85SNPfilter

# After running PLINK, rename the output file. Also, remove the rows of allele names (only because these caused issues when being read into R) and pop names.

sed 1,2d Pine_NaturalComp_85SNPfilter.recode.strct_in > Pine_NaturalComp_85SNPfilter.stru
