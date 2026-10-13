# Overview
This repository contains the code used for Goal #1 of the Morton Arboretum [IMLS 2022 National Leadership Grant (NLG)](https://www.imls.gov/grants/awarded/mg-251613-oms-22), 
which seeks to determine if, and in what circumstances, measures of geographic and ecological coverage predict genetic coverage.

To do this, we utilize datasets for which genetic and geographic data exists for **wild** individuals of different plant species, and 
we utilize resampling approaches to randomly generate subsets of individuals and measure how well those subsets reflect the total genetic, geographic, 
or ecological diversity of the complete sample set (what we term "coverage"). Further explanations of the methodology are provided below.

Several functions (specifically the underlying functions used to calculate geographic
and ecological coverage) were adapted from gap analysis approaches (for instance, see the repository [here](https://github.com/eb-bruns/conservation-gap-analysis)).

# Repository Structure
There are two primary folders in this repository: `Scripts`, which contains all the scripts used to run analyses, and `Datasets`,
which includes some of the primary and intermediate data files used to run those analyses. An explanation of each folder is provided
below.

## [`Scripts`](https://github.com/HobanLab/GeographicGeneticCorrelation/tree/main/Scripts) Folder
0. The `dataCleaning` folder contains the BASH scripts needed to process the genetic data of three of the datasets.
1. `Aim_1` contains scripts which were used for sampling strategy comparisons. The goal of these is to see if datasets designed to maximize
geographic coverage have higher genetic metric values than randomly sampled datasets. 
  - To do this, geographically maximized datasets needed to be generated for each species; this is achieved 
  in the [`runGeoSampling.R`](https://github.com/HobanLab/GeographicGeneticCorrelation/blob/main/Scripts/Aim_1/geoSamplingFunctions.R)
  script (which calls the [`geoSamplingFunctions.R`](https://github.com/HobanLab/GeographicGeneticCorrelation/blob/main/Scripts/Aim_1/geoSamplingFunctions.R) script).
  - Once geographically maximized sample sets have been generated for each species at each buffer size (CSVs in `GeoCoreSets` folders),
  these are read in for each `AltGenMetrics.R` script written for each species. These scripts are broken into various sections
    - Read in genetic/geographic data and prepare for processing
    - Calculate genetic metrics for geographically maximized datasets
      - Allelic coverage (using custom code, from [`functions_GeoGenCoverage.R`](https://github.com/HobanLab/GeographicGeneticCorrelation/blob/main/Scripts/functions_GeoGenCoverage.R))
      - Proportion of heterozygosity (HE) and average entry-to-nearest-entry distance (EN) metrics, calculated using the [CoreHunter](https://www.corehunter.org/r) R package.
    - Calculate genetic metrics for randomly sampled datasets
    - Plot differences between geographically maximized and randomly sampled approaches for each genetic metric
2. `Aim_2` contains scripts used to determine the optimal geographic (SDM approach and "Total buffer" approach) and ecological buffer sizes for each taxon.
  - This analysis for the project was originally titled "Spatial Metrics for Buffer Optimization"; hence, each script for each species is labeled "SMBO".
  - Optimal buffer sizes were determined using allelic coverage, and the functions called in these scripts are declared in the [`functions_GeoGenCoverage.R`](https://github.com/HobanLab/GeographicGeneticCorrelation/blob/main/Scripts/functions_GeoGenCoverage.R) script.
  - These scripts also contain commands for plotting the outputs of these SMBO runs. Note that none of these graphs are ultimately included in the manuscript.
3. `Aim_3` contains the [`runPointSummaries.R`](https://github.com/HobanLab/GeographicGeneticCorrelation/blob/main/Scripts/Aim_3/runPointSummaries.R) script, which was used to
build and plot the correlations between optimal buffer sizes and spatial summary metrics. 
  - This script reads in the SMBO data generated under `Aim_2`.

## [Datasets](https://github.com/HobanLab/GeographicGeneticCorrelation/tree/main/Datasets) Folder
This folder is grouped by species (which are labeled using four letters: first two letters of genus, first two letters of species). 
Within each folder are certain subfolders
1. `Genetic` includes the genetic source files for each species.
2. `Geographic` includes CSVs of occurrences of each individuals.
3. `resamplingData` includes the outputs of the analyses above: for Aim #1, `_AltGenMets.Rdata` files, and for Aims #2 and #3, SMBO resampling arrays.
  - *Note 1*: for rare species--AMTH, COGL, HIWA--only the `resamplingData` folder is present, as genetic and geographic data for these sensitive species is protected.
  - *Note 2*: some species (ARTH, PICO) have massive genetic datasets which exceed GitHub's file size limit (100MB), and so can't be linked on this repository.

# SMBO (Aim #2) Code structure
This section describes the structure of the code used to determine optimal buffer sizes for each species (Aim #2). 

## Overview
The resampling arrays (see [Outputs](https://github.com/HobanLab/GeographicGeneticCorrelation#outputs) below) are generated using a series of nested functions iterated using `sapply`. The functions at the uppermost level (`geo.gen.Resample` and its parallelized version, `geo.gen.Resample.Par`) are called in the scripts for each species. These functions are wrappers: they `sapply` the `exSituResample` function (or `exSituResample.Par`, the parallelized version) over the specified number of resampling replicates. 
In turn, `exSituResample` and `exSituResample.Par` are wrappers that use `sapply` (or `parSapply`) to reiterate the `calculateCoverage` function for every number of samples included in a wild dataset, starting at 2 and ranging all the way up to the total number of samples.

`calculcateCoverage` is the "core function" of the code structure. It is divided into sections that calculate the 
coverage values (genetic, geographic, and/or ecological) of a randomly selected subset of samples (variable name `samp`) using "worker" functions. The genetic section uses the worker function `gen.getAlleleCategories`; the geographic section uses `geo.compareBuff` and `geo.compareBuffSDM`; and the ecological section uses `eco.compareBuff`. Beyond these, there are a couple lower level functions, 
which are used for the geographic/ecological coverage calculations (`eco.intersectBuff`,`geo.checkSDMres`, and `createBuffers`).

## Inputs
The most important arguments provided to the resampling functions (`geo.gen.Resample` and `geo.gen.Resample.Par`) are:
1. a `data.frame` with 3 columns: sample name, latitude, and longitude. Lat/longs need to be in decimal degree format, and need to have the column names `decimalLatitude` and
`decimalLongitude`
2. a `genind` file, in which the order and the names of samples must match the order/names of samples in the coordinate data.frame (#1)

An error will be (intentionally) thrown if sample names/order do not match exactly between these two arguments!

Additionally, the functions require the specification of geographic/ecological buffer sizes (either a single integer or a vector of values), the Spatial Vectors representing the .shp files of polygons representing both national borders and ecoregion data (if available), and the number of resampling replicates. 

## Outputs
The uppermost resampling functions (`geo.gen.Resample` and `geo.gen.Resample.Par`) generate a single 3 dimensional array, with the dimensions as follows:
- **rows**: number of randomly selected samples, for which genetic, geographic, and ecological coverage is calculated. The first row corresponds to 2 samples, and the last row to the total number of samples.
- **columns**: coverage values of different metrics: 
	- **Column 1**: the Total allelic representation (all categories of alleles)
	- **Columns 2--5**: the representation values for alleles of different frequency categories (Very common, Common, Low frequency, Rare) 
	- **Columns 6 and up**: the geographic coverage values (for each specified buffer size) followed by the ecological coverage values (for each specified buffer size)
- **slices**: each array slice (3rd dimension) represents a different, independent resampling replicate. Averaging results across resampling replicates allows us to calculate summary statistics.

# Data and Contact
The links below direct to the repositories for the raw genetic/geographic data utilized for some of the taxa assessed. 
Certain files within each of these repositories were processed to provide the input files available in the [`Datasets`](https://github.com/HobanLab/GeographicGeneticCorrelation/tree/main/Datasets) folder.
Note that five other datasets mentioned in the manuscript came from coauthors and are not publicly available (hence, they are not linked here).

1. _Quercus acerifolia_
	+ [SRA link](https://submit.ncbi.nlm.nih.gov/subs/sra/SUB10415299/overview)
2. _Yucca brevifolia_
	+ [Dryad link](https://datadryad.org/stash/dataset/doi:10.5061/dryad.7pj4t)
3. _Quercus lobata_
	+ [Dryad link](https://datadryad.org/stash/dataset/doi:10.5061/dryad.5dv41ns4n)
4. _Pinus contorta_
	+ [Dryad link](https://datadryad.org/stash/dataset/doi:10.5061/dryad.ncjsxkstp)
5. _Mimulus guttatus_
	+ [Dryad link](https://datadryad.org/stash/dataset/doi:10.5061/dryad.ncjsxkstp)

For questions about these datasets, how they were processed, or the scripts included here, open an Issue or contact [Austin Koontz](https://akoontz11.netlify.app/).
