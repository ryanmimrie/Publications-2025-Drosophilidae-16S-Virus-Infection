<img src="img/logo.png" alt="Logo" width="200"/><br>

# *Drosophilidae* Microbiome Response to Virus Infection

This repository contains data and code used in Imrie et al., (2025) "Within and Across-Species Patterns in the *Drosophilidae* Microbiome Response to Virus Infection."

## Contents
- `data`: Contains experimental and reference dataset files.
- `figures`: Contains raw ggplot `svg` files, Affinity Designer vector files, and final `jpg` files for each manuscript figure.
- `models`: Contains `Rdata` files of high-iteration MCMCglmm runs summarised in the manuscript.
- `img`: Contains images used in this repository
- `reads`: Contains read quality profile plots (Zenodo downloads of read files should be placed here).
- `scripts`: Contains all processing, analysis, and plotting scripts used in this study.

## Software Requirements
### Platform Support
![Windows](https://img.shields.io/badge/Windows-blue?logo=microsoftwindows)
![macOS](https://img.shields.io/badge/macOS-black?logo=apple)
![Linux](https://img.shields.io/badge/Linux-grey?logo=linux)

### Dependencies
![R Version](https://img.shields.io/badge/R-4.5.0-blue)
![devtools](https://img.shields.io/badge/devtools-2.4.5-ff69b4)
![tidyverse](https://img.shields.io/badge/tidyverse-2.0.0-blue)

### System Requirements
![RAM](https://img.shields.io/badge/minimum%20RAM-8GB-important)

### Installation

To run this project, you need to have R and RStudio installed:

- [Download R](https://cran.r-project.org/)
- [Download RStudio](https://posit.co/download/rstudio-desktop/)

Within RStudio, all library dependencies can be installed by running the `scripts/setup.R` script.
