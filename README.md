# Analysis Code for Yaron et al. (2025) PLOS Biology
**"Auditory cortex neurons that encode negative prediction errors respond to omissions of sounds in a predictable sequence"**

[![MATLAB](https://img.shields.io/badge/MATLAB-R2020a%2B-orange.svg)](https://www.mathworks.com/products/matlab.html)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

This repository contains the **MATLAB analysis code** used to generate figures and results in our PLOS Biology paper. These scripts are designed specifically for our experimental dataset and are provided for **reproducibility and transparency** of our published findings.

**Associated Publication:**
> Yaron, A., Shiramatsu-Isoguchi, T., Kern, F. B., Ohki, K., Takahashi, H., & Chao, Z. C. (2025). Auditory cortex neurons that encode negative prediction errors respond to omissions of sounds in a predictable sequence. *PLOS Biology*. (Manuscript ID: PBIOLOGY-D-25-00849R2)


## Quick Start (to reproduce our results)

```matlab
% 1. Run batch processing on our 18 datasets
batchall

% 2. Identify PEONs using our criteria
FindPEONS  % Figure 2C, 2D

% 3. Generate paper figures
depthnewest    % Figure 2E-F
plot4A         % Figure 4A  
omissionvstone % Figure 4C
```

## System Requirements

### MATLAB Version
- **Minimum:** MATLAB R2020a
- **Recommended:** MATLAB R2022a or later

### Required Toolboxes
- **Statistics and Machine Learning Toolbox** (essential)
  - Functions used: `corr`, `signrank`, `fitlm`, `fitlme`, `nlinfit`



## Our Analysis Pipeline

This workflow processes our specific dataset of 18 recording sessions from rat auditory cortex:

```
Our Raw Data → Statistical Tests → Matrix Format → PEON Identification → Paper Figures
     ↓              ↓                 ↓               ↓                    ↓
  res*.mat → batchall.m → allmats.m → FindPEONS.m → [figure scripts]
  bl*.mat      ↓            ↓           ↓              ↓
tables.mat     ↓            ↓           ↓         Paper Figures:
           18 datasets  allommat    123 PEONs     2C,2D, 2E, 2F, 3A, 4A, 4C
           processed      ↓           ↓              
                     990 neurons  preferred_tone   
                     analyzed      _direction
```

## Script Reference

### Core Pipeline Scripts (run in order)

| Script | Purpose | Key Outputs | Runtime* |
|--------|---------|-------------|----------|
| `batchall.m` | Process 18 datasets, statistical tests | Workspace variables per dataset | ~60 min |
| `allmats.m` | Convert to 3D matrices (50 trials/condition) | `allommat`, `allAmat`, `allBmat` | ~2 min |
| `FindPEONS.m` | **Main analysis:** Identify PEONs using split-half | `PEONs_training`, `preferred_tone_direction` | ~5 min |

*Runtime estimates for typical dataset

### Figure Generation Scripts (run after FindPEONS.m)

| Script | Figures Generated | Description | Runtime* |
|--------|-------------------|-------------|----------|
| `depthnewest.m` | **2E, 2F** | Laminar & areal distribution | ~2 min |
| `plot4A.m` | **4A** | Population responses vs probability | ~3 min |
| `omissionvstone.m` | **4C** | PEON selectivity pie charts | ~2 min |
| `rampuphalf2.m` | **3A** | Trial-by-trial buildup dynamics | ~3 min |

*Runtime estimates for typical dataset

### Support Scripts

| Script | Purpose |
|--------|---------|
| `maketests1214NosideAllAud.m` | Process single dataset (called by `batchall.m`) |
| `createtable.m` | Aggregate data across datasets |

## Key Parameters

### FindPEONS.m Configuration
```matlab
ANALYSIS_MODE = 'ODD';     % Options: 'ODD', 'EVEN', 'ALL'
ALPHA_CORR = 0.05;         % Significance level for correlations  
ALPHA_WILCOX = 0.05;       % Significance level for Wilcoxon tests
```

### Our Specific Dataset
- **Recording sessions:** 18 penetrations in rat auditory cortex
- **Total neurons:** ~990 analyzed
- **PEONs identified:** 123 neurons  
- **Trials per condition:** 50 (after subsampling)
- **Probability conditions:** 8 levels (0%, 5%, 10%, 20%, 75%, 85%, 90%, 95%)
- **Brain areas:** A1, VAF, AAF

## Outputs 

### Paper Figures Generated
- **Figure 2E, 2F:** Laminar and areal distribution of PEONs (`depthnewest.m`)
- **Figure 3A:** Trial-by-trial response buildup (`rampuphalf2.m`) 
- **Figure 4A:** Population responses vs probability (`plot4A.m`)
- **Figure 4C:** PEON selectivity analysis (`omissionvstone.m`)

### Generated Files
- **PDF files:** High-resolution figures matching those in paper
- **MAT files:** Intermediate analysis results and workspace variables

### Key Variables
```matlab
PEONs_training          % Indices of identified PEONs
preferred_tone_direction % Direction preference (+1/-1) per PEON
allommat               % 3D matrix: neurons × trials × conditions
T                      % Table with neuron metadata (depth, location)
```

## Data Structure

### Input Files
```
tables.mat              % Dataset metadata (18 datasets)
res[ID].mat            % Raw neural responses  
bl[ID].mat             % Baseline activity
```


## Code Availability
This code is archived on Zenodo 
- **GitHub (development):** https://github.com/Amityar01/PEON/ 
- **Zenodo (permanent archive):** 

## Data Availability
- **Raw experimental data:** Available on Zenodo https://doi.org/10.5281/zenodo.15531781


If you're working with different auditory cortex data, you would need to substantially modify these scripts or develop your own analysis approach.


## License

This project is licensed under the MIT License - see [LICENSE](LICENSE) file for details.

## Support & Contact

**For reproducing our paper results:**
- **Issues with our code:** Use GitHub Issues  
- **Questions about our methods:** Contact Amit Yaron at amityar@gmail.com
- **Data access problems:** Check Zenodo DOI above


## Contributors

- **Amit Yaron** - Primary developer and analysis design
- **Zenas C. Chao** - Senior supervision and methodology

---
