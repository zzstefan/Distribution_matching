# Distribution Matching for Brain Connectivity Harmonization

## Overview
This repository contains the implementation of a distribution-matching approach to harmonizing structural brain connectivity across different sites and scanners. The method aims to align the statistical properties of connectivity data from different datasets, improving the reliability and comparability of structural connectivity studies that combine data from multiple sources.

## Background
Multi-site diffusion-weighted magnetic resonance imaging (dMRI) studies offer enhanced statistical power for investigating brain structure, but face challenges due to variations in scanner hardware and acquisition protocols. This distribution-matching technique addresses the harmonization of structural brain connectivity data by aligning the statistical properties of connectivity from different datasets.

## Features
- Gamma distribution modeling of structural connectivity values
- Distribution matching to align data from different scanners/sites
- Comparison with ComBat/Covbat harmonization method
- Evaluation through correlation analysis with clinical scores (MMSE) and demographic variables (age)

## Installation

### Prerequisites
- MATLAB (tested with R2020b and later)
- FreeSurfer (for processing anatomical MR images)
- FSL (FMRIB Software Library)
- CSAODF (Diffusion MRI Orientation Distribution Function in Constant Solid Angel and Hough-Transform Tractography, https://www.nitrc.org/projects/csaodf-hough)

### Setup
1. Clone this repository:

## Usage

### Data Preparation
The code is designed to work with connectivity matrices generated from dMRI data. For data preprocessing:
1. Process raw dMRI data using FreeSurfer and FSL
2. Generate structural connectivity matrices using CSAODF package

### Running the Distribution Matching

1. Two-Site Harmonization
    In the '/Distribution matching/Two sites/' folder, experiments focus on harmonizing data from ADNI-2 to match OASIS-3 (reference site). This includes:

   - Standard implementation (combining all subjects)
   - Sex-separated implementation (performing harmonization separately for males and females)

2. Three-Site Harmonization
    In the 'Distribution matching/Three sites/' folder, experiments include multiple reference sites:

    - Using OASIS-3 as the reference site
    - Using PreventAD as the reference site
    - Both with standard and sex-separated implementations