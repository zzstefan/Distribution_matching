# Distribution Matching for Brain Connectivity Harmonization

## Overview
This repository contains the implementation of a distribution-matching approach to harmonizing structural brain connectivity across different sites and scanners. The method aims to align the statistical properties of connectivity data from different datasets, improving the reliability and comparability of structural connectivity studies that combine data from multiple sources.

## Background
Multi-site diffusion-weighted magnetic resonance imaging (dMRI) studies offer enhanced statistical power for investigating brain structure, but face challenges due to variations in scanner hardware and acquisition protocols. This distribution-matching technique addresses the harmonization of structural brain connectivity data by aligning the statistical properties of connectivity from different datasets.

## Features
- Gamma distribution modeling of structural connectivity values
- Distribution matching to align data from different scanners/sites
- Comparison with ComBat harmonization method
- Evaluation through correlation analysis with clinical scores (MMSE) and demographic variables (age)

## Installation

### Prerequisites
- MATLAB (tested with R2020b and later)
- FreeSurfer (for processing anatomical MR images)
- FSL (FMRIB Software Library)

### Setup
1. Clone this repository: