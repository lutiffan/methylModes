[![DOI](https://zenodo.org/badge/680979978.svg)](https://doi.org/10.5281/zenodo.17448517)

# methylModes

A computationally efficient R package for detecting multimodal distributions in DNA methylation data,  implementing a peak detection algorithm based on kernel density estimation. MethylModes wraps a Shiny application for interactive analysis and includes batch processing functionality for genome-scale epigenetic studies. 

## Overview

MethylModes analyzes DNA methylation array data. It identifies probes for which the sample distributions follow multimodal distributions. This can be caused by biological variation such as probe sequence overlaps with genetic variants or technical artifacts like probe cross-reactivity.

## Installation

### Prerequisites

- R version 4.2.2 or higher
- RStudio (recommended for interactive use)

### Install methylModes

#### Method 1: using devtools
```r
# First ensure that devtools is installed
if (!require("devtools")) install.packages("devtools")
```
Make sure to run the above code by itself, and that the installation is successful before trying to execute any other lines of code.
```r
# Install from GitHub
devtools::install_github("lutiffan/methylModes")

# Verify installation by checking file location
system.file(package="methylModes")
```
#### Method 2 (Linux/Mac only): run bash script

1. Download the script (install_methylModes.sh) directly from this Github repo
2. Open a terminal window and navigate to the directory where you installed the script
3. Run `bash install_methylModes.sh`
4. The script will alert you of missing requirements, including minimum R version and compilers.

#### Troubleshooting Suggestions

* We have found during testing that when using `devtools::install_github()`, when prompted to update package dependencies, skipping updates can result in the installation becoming stuck indefinitely. Selecting option "1: All" and allowing all packages to update may resolve the issue.

* Windows: unlike Mac and Linux users, Windows users need to ensure that RTools is installed so they can compile C, C++, or Fortran source code.
Check for RTools by running
```r
if (!require("pkgbuild")) install.packages("pkgbuild")
find_rtools()
```
If you don't have it, install it from here: [Toolchains for building R and R packages from source on Windows](https://cran.r-project.org/bin/windows/Rtools/)

* Mac: ensure that XCode command line tools is installed. Open Terminal and run `xcode-select --version`.
If no version is found, [here are instructions to install it.](https://www.freecodecamp.org/news/install-xcode-command-line-tools/)

* If, after installing the package, you see an error like "Error in fetch(key) : lazy-load database ... methylModes/help/methylModes.rdb' is corrupt" then try restarting your R session and reloading the package.

## Running the Shiny Application

### Method 1: Using the Package Function

```r
library(methylModes)

# Launch the Shiny app (opens in default browser)
runMethylModesApp()

```

### Method 2: Direct Shiny App Launch

```r
# Navigate to the package directory and run
shiny::runApp("inst/shiny")

```

## Using the Test Dataset

The package includes a toy dataset for testing and demonstration purposes:

```r
# Load the test dataset
data(toyMultimodalData)

# View the dataset structure
str(toyMultimodalData)
dim(toyMultimodalData)  # 10000 probes × 800 samples

# Access the data
head(toyMultimodalData)
```

### Toy Dataset Description

- **Size**: 10,000 probes × 800 samples
- **Format**: Matrix with values ranging from 0 to 1, drawn from a Gaussian mixture distribution
- **Probe IDs**: First 10000 from IlluminaManifest450k IlmnID column ("cg00035864" ... "cg23571457")
- **Sample IDs**: "subject1" through "subject800"
- **Purpose**: Verify Shiny app functionality with multimodal distributions

### Using Test Data in the Shiny App

0. Find package file locations by running `system.file(package="methylModes")`
1. Launch the Shiny app using `runMethylModesApp()`
2. In the "Get Started" section, locate the included toy data by clicking "Browse," navigating to the package directory, and opening the "data" directory
3. Select toyMultimodalData.rda
4. You can then run analyses on individual probes or perform batch processing

## Data Format Requirements

### Input Data Format

Your methylation data should be in one of these formats:

- **CSV file**: Comma-separated values with probe IDs as row names
- **RDS file**: R data object (matrix or data.frame)
- **RDA/RData file**: R data archive
- **QS file**: Quick serialization format

### Data Structure

Inputs should be a matrix with named rows, or a data structure that can be coerced to matrix.

- **Rows**: Rownames are required. CpG probe IDs (e.g., "cg00000029")
- **Columns**: Sample IDs
- **Values**: Beta values between 0 and 1
- **Missing values**: NA or empty cells are acceptable

### Example Data Structure

```
          Sample1  Sample2  Sample3
cg00000029    0.45    0.52    0.48
cg00000108    0.78    0.81    0.79
cg00000109    0.23    0.19    0.25
```

## Methylation Array Annotation Datasets

The package includes several annotation datasets:

- `IlluminaManifest450k`: 450k array annotation data
- `IlluminaManifestEPIC`: EPIC v1.0 array annotation data  
- `IlluminaManifestEPICv2`: EPIC v2.0 array annotation data

## Features

### General
- Interactive data upload and validation
- On-demand individual probe analysis with customizable parameters
- Visualization allows preview or post-hoc analysis of how hyperparameter choice affects analysis outcomes
- Batch processing for multiple probes, with the option to subset by genomic region
- Parallel processing for multiple-probe analysis
- Results export in CSV format
- Downloaded results can be re-uploaded to the app to browse results without re-running analysis

## Citation

If you use methylModes in your research, please cite us using:

```r
citation("methylModes")
```

## Support

For issues, questions, or feature requests:
- Email: lutiffan@umich.edu
- GitHub Issues: [Create an issue](https://github.com/lutiffan/methylModes/issues)

## License

This package is licensed under GPL (>= 3). See the LICENSE file for details. 
