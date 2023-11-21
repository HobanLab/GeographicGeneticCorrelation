# Overview
This repository contains the code used for Goal #1 of the Hoban [IMLS 2022 National Leadership Grant (NLG)](https://www.imls.gov/grants/awarded/mg-251613-oms-22), 
which seeks to determine if, and in what circumstances, measures of geographic and ecological coverage _ex situ_ are predictive of genetic coverage.

To do this, we utilize datasets for which genetic and geographic data exists for individuals (or populations), and 
we utilize resampling approaches to randomly generate subsets of individuals and measure how well those subsets reflect the total genetic, geographic, 
or ecological diversity of the complete sample set (what we term "coverage"). We iterate this process for subsets of increasing sizes, and repeat this resampling approach
multiple times to account for the stochasticity of random sampleing. Summary statistics of coverage metrics across resampling replicates are then calculated.

The methodology for this project is very much in development. Several functions (specifically the underlying functions used to calculate geographic
and ecological coverage) were adapted from gap analysis approaches (for instance, see the repository [here](https://github.com/eb-bruns/conservation-gap-analysis)).

# Repository layout
The functions used for this project build off of one another and are inteded to provide a simplified interface at the "uppermost" level of the code 
(i.e. the level at which data is read in). All functions are declared in the [`functions_GeoGenCoverage.R`](https://github.com/HobanLab/GeographicGeneticCorrelation/blob/main/functions_GeoGenCoverage.R)
script. 

There's currently a single analysis script for each species analyzed, but this layout might change as the project develops and more datasets are processed.

## Code structure
Resampling arrays (see [Outputs](https://github.com/HobanLab/GeographicGeneticCorrelation#outputs) below) are generated using a series of nested functions 
iterated using `sapply`. The functions at the uppermost level (`geo.gen.Resample` and its parallelized version, `geo.gen.Resample.Parallel`) are called in 
the scripts for each species. These functions are wrappers: they `sapply` the `exSituResample` function over the specified number of resampling replicates. 
In turn, `exSituResample` is a wrapper that uses `sapply` to reiterate the `calculateCoverage` function for every number of samples included in a wild dataset, 
starting at 2 and ranging all the way up to the total number of samples.

`calculcateCoverage` is the "core function" of the code structure. It is divided into sections that calculate the 
coverage values (genetic, geographic, and/or ecological) of a subset of randomly selected samples (variable name `samp`) using "worker" functions. 
The genetic section uses the worker function `gen.getAlleleCategories`; the geographic section uses `geo.compareBuff`; and the ecological section uses `eco.compareBuff`. Beyond these, there are a couple lower level functions, 
which are used for the geographic/ecological coverage calculations (`eco.intersectBuff` and `createBuffers`).

## Inputs
The most important arguments provided to the resampling functions (`geo.gen.Resample` and `geo.gen.Resample.Parallel`) are:
1. a `data.frame` with 3 columns: sample name, latitude, and longitude. Lat/longs need to be in decimal degree format, and need to have the column names `decimalLatitude` and
`decimalLongitude`
2. a `genind` file, in which the order and the names of samples must match the order/names of samples in the coordinate data.frame (#1)

An error will be (intentionally) thrown if sample names/order do not match exactly between these two arguments!

Additionally, the functions require the specification of geographic/ecological buffer sizes, the Spatial Vectors representing the .shp files of polygons
representing both national borders and ecoregion data (if available), and the number of resampling replicates. 

## Outputs
The uppermost resampling functions (`geo.gen.Resample` and `geo.gen.Resample.Parallel`) generate a single 3 dimensional array, with the dimensions as follows:
- **rows**: number of randomly selected samples, for which genetic, geographic, and ecological coverage is calculated. The first row corresponds to 2 samples, and the last row to the total number of samples.
- **columns**: coverage values of different metrics: 
	- **Column 1**: the Total allelic representation (all categories of alleles)
	- **Columns 2--5**: the representation values for alleles of different frequency categories (Very common, Common, Low frequency, Rare) 
	- **Column 6**: the geographic coverage values (optional) 
	- **Column 7**: the ecological coverage values (optional)
- **slices**: each array slice (3rd dimension) represents a different, independent resampling replicate. Averaging results across resampling replicates allows us to calculate summary statistics.

## Analysis
After building resampling arrays, we pass the results into linear models which specify allelic representation values as the response variable, and either geographic
or ecological coverage values as the explanatory variable. We capture the R-squared values generated from linear models, and create 2 different types of plots to illustrate
the relationships between the coverage estimates:
- "correlation plots", which plot the average genetic coverage values (y-axis) versus the average geographic/ecological coverage values (x-axis)
- "coverage plots", where the average coverage values for all 3 measures are plotted in different colors against the number of samples in the dataset

These plots are generated in the 2nd half of the demo scripts, for each species.

# Data and Contact
## Datasets
1. _Quercus acerifolia_
	+ [SRA link](https://submit.ncbi.nlm.nih.gov/subs/sra/SUB10415299/overview)
2. _Quercus lobata_
	+ [Dryad link](https://datadryad.org/stash/dataset/doi:10.5061/dryad.5dv41ns4n)
3. _Pinus contorta_
	+ [Dryad link](https://datadryad.org/stash/dataset/doi:10.5061/dryad.ncjsxkstp)
4. _Mimulus guttatus_
	+ [Dryad link](https://datadryad.org/stash/dataset/doi:10.5061/dryad.ncjsxkstp)


For questions about this datasets or the scripts included here, open an Issue or contact [Austin Koontz](https://akoontz11.netlify.app/).
