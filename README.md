# Distribution Matching for Brain Connectivity Harmonization

## Overview
This repository contains the implementation of a distribution-matching approach to harmonizing structural brain connectivity across different sites and scanners. The method aims to align the statistical properties of connectivity data from different datasets, improving the reliability and comparability of structural connectivity studies that combine data from multiple sources.

## Background
Multi-site diffusion-weighted magnetic resonance imaging (dMRI) studies offer enhanced statistical power for investigating brain structure, but face challenges due to variations in scanner hardware and acquisition protocols. This distribution-matching technique addresses the harmonization of structural brain connectivity data by aligning the statistical properties of connectivity from different datasets.

## Features
- Gamma+delta distribution modeling of structural connectivity values
- Distribution matching to align data from different scanners/sites
- Comparison with ComBat/CovBat harmonization method
- Evaluation through correlation analysis with clinical scores (MMSE) and demographic variables (age)

## Installation

### Optional Prerequisites
- ComBat (https://github.com/rpomponio/neuroHarmonize)
- CSAODF (Diffusion MRI Orientation Distribution Function in Constant Solid Angel and Hough-Transform Tractography, https://www.nitrc.org/projects/csaodf-hough)

### Setup
1. Clone this repository.

## Usage

### Data Preparation
The code is designed to work with connectivity matrices generated from dMRI data. For data preprocessing:
1. Process raw dMRI data, e.g. using FreeSurfer and FSL
2. Generate structural connectivity matrices, e.g. using CSAODF package

### Running the Distribution Matching

#### Two-Site Harmonization (OASIS-3 and ADNI-2)
In the `/Distribution matching/Two sites/` folder, experiments focus on harmonizing data from ADNI-2 to match OASIS-3 (reference site). This includes:

- Standard implementation (combining all subjects)
- Sex-separated implementation (performing harmonization separately for males and females)

#### Three-Site Harmonization (OASIS-3, ADNI-2 and PREVENT-AD)
In the `/Distribution matching/Three sites/` folder, experiments include multiple reference sites:

- Using OASIS-3 as the reference site
- Using PREVENT-AD as the reference site
- Both with standard and sex-separated implementations

### Correlation analysis
Our method incorporates comprehensive correlation analysis to evaluate harmonization quality:

#### MMSE Score Correlation
In the `/Correlation analysis/` folder, `run_MMSE_correlation.py`
  - Measures correlation between MMSE and structural connectivity before and after harmonization

#### Age Correlation
In the `/Correlation analysis/` folder, `run_age_correlation.py`
  - Analyzes correlation between age and structural connectivity before and after harmonization

## Datasets

The method was validated using:

- **OASIS-3** (Open Access Series of Imaging Studies)
- **ADNI-2** (Alzheimer's Disease Neuroimaging Initiative)
- **PREVENT-AD** (PResymptomatic EValuation of Experimental or Novel Treatments for AD)

*Note: The datasets are not included in this repository due to their size and licensing restrictions. They can be obtained from their respective official sources.*

## Citation
If you use this code in your research, please cite our paper:
```
Zhou Z, Fischl B, Aganj I, Alzheimer's Disease Neuroimaging Initiative. Harmonization of Structural Brain Connectivity Through Distribution Matching. Hum Brain Mapp. 2025 Jun 15; 46(9):e70257. PMID: 40545448
```

## Acknowledgments

This work was supported by Michael J. Fox Foundation for Parkinson‘s Research, MJFF- 021226. Additional funding was provided by the National Institute on Aging, RF1AG068261, R01AG068261, R01AG008122, R01AG016495, R01AG064027, R01AG070988, R21AG082082 National Institute of Biomedical Imaging and Bioengineering, P41EB030006, R01EB006758, R01EB019956, R01EB023281, R21EB018907 the Brain Initiative Brain Connects consortium, U01NS132181, UM1NS132358 National Institute of Neurological Disorders and Stroke, R01NS0525851, R01NS070963, R01NS083534, R01NS105820, R21NS072652, U01NS086625, U24NS10059103 National Institute of Mental Health, U01MH117023, R01MH121885, R01MH123195, RF1MH123195, NIH Blueprint for Neuroscience Research, U01MH093765.

## Contact

For questions, issues, or collaboration opportunities, please contact:
- Zhen Zhou (zzhou22@mgh.harvard.edu)
- Iman Aganj (iaganj@mgh.harvard.edu)
