# Overview
This repository contains the code used for Goal #1 of the Morton Arboretum [IMLS 2022 National Leadership Grant (NLG)](https://www.imls.gov/grants/awarded/mg-251613-oms-22), 
which seeks to determine if, and in what circumstances, measures of geographic and ecological coverage predict genetic coverage.

To do this, we utilize datasets for which genetic and geographic data exists for **wild** individuals (or populations) of different plant species, and 
we utilize resampling approaches to randomly generate subsets of individuals and measure how well those subsets reflect the total genetic, geographic, 
or ecological diversity of the complete sample set (what we term "coverage"). We repeat this process for subsets of increasing sizes, and also generate multiple resampling replicates
to account for the stochasticity of random sampling. Summary statistics of coverage metrics across resampling replicates are then calculated.

The methodology for this project is very much in development. Several functions (specifically the underlying functions used to calculate geographic
and ecological coverage) were adapted from gap analysis approaches (for instance, see the repository [here](https://github.com/eb-bruns/conservation-gap-analysis)).

## Repository layout
There are two primary folders in this repository.

1. [`Scripts`](https://github.com/HobanLab/GeographicGeneticCorrelation/tree/main/Scripts) 
contains all of the R and BASH scripts used to analyze datasets. Included in this folder is the [`functions_GeoGenCoverage.R`](https://github.com/HobanLab/GeographicGeneticCorrelation/blob/main/Scripts/functions_GeoGenCoverage.R) script, which declares all of the R functions used in all the other scripts. The functions used for these 
analyses build off of one another and are intended to provide a simplified interface at the "uppermost" level of the code 
(i.e. the level at which data is read in). Additionally, there's a single analysis R script for each species analyzed. This layout may change as the project develops and more datasets are processed. Any other required scripts are also included in this folder.

2. [`Datasets`](https://github.com/HobanLab/GeographicGeneticCorrelation/tree/main/Datasets) contains the genetic and geographic
input files required for each species. The subfolders within this folder are organized by species: for each species, there's
a `Genetic` folder, a `Geographic` folder, and a `resamplingData` folder, where resampling arrays are saved. Note that in 
some instances, the input files used for analyses cannot be included on the GitHub repo, due to surpassing file size limits
(50 MB); in these cases, the files are included in the [`.gitignore`](https://github.com/HobanLab/GeographicGeneticCorrelation/blob/main/.gitignore) file, 
to avoid being tracked on GitHub. Shapefiles for geographic (global country borders) and ecological 
(ecoregion data layers) analyses are also not tracked in the repository, because of file size limitations.

## Code structure
Resampling arrays (see [Outputs](https://github.com/HobanLab/GeographicGeneticCorrelation#outputs) below) are generated using a series of nested functions iterated using `sapply`. The functions at the uppermost level (`geo.gen.Resample` and its parallelized version, `geo.gen.Resample.Par`) are called in the scripts for each species. These functions are wrappers: they `sapply` the `exSituResample` function (or `exSituResample.Par`, the parallelized version) over the specified number of resampling replicates. 
In turn, `exSituResample` and `exSituResample.Par` are wrappers that use `sapply` (or `parSapply`) to reiterate the `calculateCoverage` function for every number of samples included in a wild dataset, starting at 2 and ranging all the way up to the total number of samples.

`calculcateCoverage` is the "core function" of the code structure. It is divided into sections that calculate the 
coverage values (genetic, geographic, and/or ecological) of a randomly selected subset of samples (variable name `samp`) using "worker" functions. The genetic section uses the worker function `gen.getAlleleCategories`; the geographic section uses `geo.compareBuff`; and the ecological section uses `eco.compareBuff`. Beyond these, there are a couple lower level functions, 
which are used for the geographic/ecological coverage calculations (`eco.intersectBuff` and `createBuffers`).

### Inputs
The most important arguments provided to the resampling functions (`geo.gen.Resample` and `geo.gen.Resample.Par`) are:
1. a `data.frame` with 3 columns: sample name, latitude, and longitude. Lat/longs need to be in decimal degree format, and need to have the column names `decimalLatitude` and
`decimalLongitude`
2. a `genind` file, in which the order and the names of samples must match the order/names of samples in the coordinate data.frame (#1)

An error will be (intentionally) thrown if sample names/order do not match exactly between these two arguments!

Additionally, the functions require the specification of geographic/ecological buffer sizes, the Spatial Vectors representing the .shp files of polygons representing both national borders and ecoregion data (if available), and the number of resampling replicates. 

### Outputs
The uppermost resampling functions (`geo.gen.Resample` and `geo.gen.Resample.Par`) generate a single 3 dimensional array, with the dimensions as follows:
- **rows**: number of randomly selected samples, for which genetic, geographic, and ecological coverage is calculated. The first row corresponds to 2 samples, and the last row to the total number of samples.
- **columns**: coverage values of different metrics: 
	- **Column 1**: the Total allelic representation (all categories of alleles)
	- **Columns 2--5**: the representation values for alleles of different frequency categories (Very common, Common, Low frequency, Rare) 
	- **Column 6**: the geographic coverage values (optional) 
	- **Column 7**: the ecological coverage values (optional)
- **slices**: each array slice (3rd dimension) represents a different, independent resampling replicate. Averaging results across resampling replicates allows us to calculate summary statistics.

## Analysis
After building resampling arrays, we pass the results into linear models which specify allelic representation values as the response variable, and either geographic or ecological coverage values as the explanatory variable. We capture the R-squared values generated from linear models, and create 2 different types of plots to illustrate the relationships between the coverage estimates:
- "correlation plots", which plot the average genetic coverage values (y-axis) versus the average geographic/ecological coverage values (x-axis)
- "coverage plots", where the average coverage values for all 3 measures are plotted in different colors against the number of samples in the dataset

These plots are generated in the 2nd half of the demo scripts, for each species.

# Data and Contact
The links below direct to the repositories for the raw genetic/geographic data utilized for each species. Certain files within each of these repositories were processed to provide the input files available in the [`Datasets`](https://github.com/HobanLab/GeographicGeneticCorrelation/tree/main/Datasets) folder.

1. _Quercus acerifolia_
	+ [SRA link](https://submit.ncbi.nlm.nih.gov/subs/sra/SUB10415299/overview)
2. _Quercus lobata_
	+ [Dryad link](https://datadryad.org/stash/dataset/doi:10.5061/dryad.5dv41ns4n)
3. _Pinus contorta_
	+ [Dryad link](https://datadryad.org/stash/dataset/doi:10.5061/dryad.ncjsxkstp)
4. _Mimulus guttatus_
	+ [Dryad link](https://datadryad.org/stash/dataset/doi:10.5061/dryad.ncjsxkstp)
5. _Yucca brevifolia_
	+ [Dryad link](https://datadryad.org/stash/dataset/doi:10.5061/dryad.7pj4t)

For questions about these datasets, how they were processed, or the scripts included here, open an Issue or contact [Austin Koontz](https://akoontz11.netlify.app/).
