# Overview
This repository contains the code used for Goal #1 of the Hoban [IMLS 2022 National Leadership Grant (NLG)](https://www.imls.gov/grants/awarded/mg-251613-oms-22), 
which seeks to determine when, and in what circumstances, measures of geographic and ecological coverage _ex situ_ can be used as proxies for levels of genetic coverage.

To answer this question, we utilize datasets for which there genetic, geographic, and ecological data for individuals (or populations), and 
we utilize resampling approaches to randomly generate subsets of individuals and measure how well those subsets reflect the total genetic, geographic, 
or ecological diversity of the complete sample set (what we term "coverage"). We iterate this process for different subset sizes to account for the 
stochasticity of random sampling, and generate summary metrics of our coverage results across replicates.

The methodology for this approach is very much in development. Several functions (specifically the underlying functions used to calculate geographic
and ecological coverage) were adapted from gap analysis approaches (for instance, see the repository [here](https://github.com/eb-bruns/conservation-gap-analysis)).

# Repository layout
The [functions](https://github.com/HobanLab/GeographicGeneticCorrelation/blob/main/functions_GeoGenCoverage.R) used for this project build off of one another and 
are inteded to provide a simplified interface at the "uppermost" level of the code (i.e. the level at which data is read in). Currently, analysis
scripts are organized by species, but this layout might change as the project develops and more datasets are analyzed.

## Code structure
Resampling arrays (see [Outputs](https://github.com/HobanLab/GeographicGeneticCorrelation#outputs) below) are generated using a series of nested functions 
iterated using `sapply`. The functions at the uppermost level (`geo.gen.Resample` and its parallelized version, `geo.gen.Resample.Parallel`) are called in 
the demo scripts for each species. These functions are wrappers: they `sapply` the `exSituResample` function over the number of specified resampling replicats. 
In turn, `exSituResample` is a wrapper that uses `sapply` to reiterate the `calculateCoverage` function for every number of samples included in a wild dataset, 
starting at 2 and ranging all the way up to the total number of samples.

`calculcateCoverage` is a wrapper of several different functions, and is the "core function" of the code structure. It is divided into sections that calculate the 
coverage values of a subset of randomly selected samples (variable name `samp`) using worker functions. The genetic section uses the worker function `getAlleleCategories`; 
the geographic section uses `geo.compareBuff`; and the ecological section uses `eco.compareBuff`. Beyond these, there are a couple lower level functions, 
used for the geographic/ecological coverage calculations (`eco.intersectBuff` and `createBuffers`).

## Inputs
The most important arguments provided to resampling functions are
1. a `data.frame` with 3 columns: sample name, latitude, and longitude. Lat/longs need to be in decimal degree format
2. a genind file in which the order and the names of samples match the order/names of samples in the coordinate data.frame (#1)

An error will be thrown if sample names/order do not match between these two arguments!

Additionally, the functions require the specification of geographic/ecological buffer sizes, the Spatial Vectors representing the .shp files of polygons
representing both national borders and ecoregion data (if available), and the number of resampling replicates. 

## Outputs
The uppermost resampling functions (`geo.gen.Resample` and its parallelized version, `geo.gen.Resample.Parallel`) generate a single 3 dimensional array,
with the dimensions as follows:
- rows: number of randomly selected samples, for which genetic, geographic, and ecological coverage is calculated
- columns: coverage values of different metrics. Column 1 is the Total allelic representation; Columns 2--5 are the allelic representation values for different
frequency categories of alleles (Very common, Common, Low frequency, Rare); Column 6 contains the geographic coverage values; and Column 7 contains the ecological 
coverage values (assuming ecological coverage calculation is specified in the resampling function call--this is optional)
- slices: each array slice represents a different resampling replicate

# Data and Contact
## Datasets
1. _Quercus acerifolia_
	+ [SRA link](https://submit.ncbi.nlm.nih.gov/subs/sra/SUB10415299/overview)
2. _Quercus lobata_
	+ [Dryad link](https://datadryad.org/stash/dataset/doi:10.5061/dryad.5dv41ns4n)
3. _Pinus contorta_
	+ [Dryad link](https://datadryad.org/stash/dataset/doi:10.5061/dryad.ncjsxkstp)
