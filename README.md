# Overview
This repository contains the code used for Goal #1 of the Hoban [IMLS 2022 National Leadership Grant (NLG)](https://www.imls.gov/grants/awarded/mg-251613-oms-22), 
which seeks to determine when, and in what circumstances, measures of geographic and ecological coverage _ex situ_ can be used as proxies for levels of genetic coverage.

To answer this question, we utilize datasets for which there genetic, geographic, and ecological data for individuals (or populations), and 
we utilize resampling approaches to randomly generate subsets of individuals and measure how well those subsets reflect the total genetic, geographic, 
or ecological diversity of the complete sample set (what we term "coverage"). We iterate this process for different subset sizes to account for the 
stochasticity of random sampling, and generate summary metrics of our coverage results.

The methodology for this approach is very much in development. Several functions (specifically the underlying functions used to calculate geographic
and ecological coverage) were adapted from gap analysis approaches (for instance, see the repository [here](https://github.com/eb-bruns/conservation-gap-analysis)).

# Repository layout
The [functions](https://github.com/HobanLab/GeographicGeneticCorrelation/blob/main/functions_GeoGenCoverage.R) used for this project build off of one another and 
are inteded to provide a simplified interface at the "uppermost" level of the code (i.e. at the level of code at which data is read in). Currently, additional 
scripts are organized by species, but this layout might change as the project develops and more datasets are analyzed.

# Data and Contact
## Datasets
1. _Quercus acerifolia_
	+ [SRA link](https://submit.ncbi.nlm.nih.gov/subs/sra/SUB10415299/overview)
2. _Quercus lobata_
	+ [Dryad link](https://datadryad.org/stash/dataset/doi:10.5061/dryad.5dv41ns4n)
