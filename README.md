# Analysis Code for Yaron et al. (2025) PLOS Biology
**"Auditory cortex neurons that encode negative prediction errors respond to omissions of sounds in a predictable sequence"**

[![MATLAB](https://img.shields.io/badge/MATLAB-R2020a%2B-orange.svg)](https://www.mathworks.com/products/matlab.html)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

This repository contains the **exact MATLAB analysis code** used to generate all figures and results in our PLOS Biology paper. These scripts are designed specifically for our experimental dataset and are provided for **reproducibility and transparency** of our published findings.

**Associated Publication:**
> Yaron, A., Shiramatsu-Isoguchi, T., Kern, F. B., Ohki, K., Takahashi, H., & Chao, Z. C. (2025). Auditory cortex neurons that encode negative prediction errors respond to omissions of sounds in a predictable sequence. *PLOS Biology*. (Manuscript ID: PBIOLOGY-D-25-00849R2)

⚠️ **Important:** These scripts are tailored to our specific experimental paradigm, data structure, and analysis requirements. They are not designed as general-purpose tools for other datasets.

## Purpose

- **Reproduce** all figures and analyses from our paper
- **Verify** our computational methods and statistical tests  
- **Provide transparency** in our analysis workflow
- **Enable replication** of our exact results

## Quick Start (to reproduce our results)

```matlab
% 1. Run batch processing on our 18 datasets
batchall

% 2. Prepare 3D matrices (50 trials/condition)
allmats

% 3. Identify PEONs using our criteria
FindPEONS

% 4. Generate paper figures
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

### Hardware
- **RAM:** 8GB minimum, 16GB recommended
- **Storage:** ~2GB for data and outputs

## Reproducing Our Results

1. **Clone this repository:**
   ```bash
   git clone https://github.com/username/yaron-peon-2025.git
   cd yaron-peon-2025
   ```

2. **Download our experimental data:**
   - Get `tables.mat` and our raw data files (`res*.mat`, `bl*.mat`) from Zenodo
   - **Data DOI:** `[TODO: Insert Zenodo DOI]`
   - Place files in the same directory as scripts

3. **Verify MATLAB setup:**
   ```matlab
   % Check required toolbox
   ver('stats')  % Should show Statistics and Machine Learning Toolbox
   ```

## Our Analysis Pipeline

This workflow processes our specific dataset of 18 recording sessions from rat auditory cortex:

```
Our Raw Data → Statistical Tests → Matrix Format → PEON Identification → Paper Figures
     ↓              ↓                 ↓               ↓                    ↓
  res*.mat → batchall.m → allmats.m → FindPEONS.m → [figure scripts]
  bl*.mat      ↓            ↓           ↓              ↓
tables.mat     ↓            ↓           ↓         Paper Figures:
           18 datasets  allommat    123 PEONs     2E, 2F, 3A, 4A, 4C
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

## Outputs (Reproduce Paper Results)

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

## Troubleshooting

### Common Issues

**"File not found" errors:**
```matlab
% Check file paths
which tables.mat
dir('res*.mat')
```

**"Undefined function" errors:**
```matlab
% Verify Statistics Toolbox
license('test', 'statistics_toolbox')
```

**Out of memory errors:**
```matlab
% Increase memory or process subsets
clear variables  % Clear workspace
memory          % Check available memory
```

### Debugging Tips
- Run scripts section by section using `%%` breaks
- Use `whos` to check workspace variables
- Enable debugging: `dbstop if error`

## Data Structure

### Input Files
```
tables.mat              % Dataset metadata (18 datasets)
res[ID].mat            % Raw neural responses  
bl[ID].mat             % Baseline activity
```

### Directory Structure
```
yaron-peon-2025/
├── README.md
├── LICENSE
├── *.m                # Analysis scripts
├── data/              # Raw data files (download separately)
│   ├── tables.mat
│   ├── resC1.mat, resC2.mat, ...
│   └── blC1.mat, blC2.mat, ...
└── outputs/           # Generated figures
    └── *.pdf          # Figure files
```

## Code Availability
This code is archived on Zenodo with DOI: **[TODO: Insert Zenodo DOI after archiving]**
- **GitHub (development):** https://github.com/username/yaron-peon-2025  
- **Zenodo (permanent archive):** [DOI to be added]

## Data Availability
- **Raw experimental data:** Available on Zenodo [DOI: TBD]
- **Analysis code:** Archived on Zenodo [DOI: TBD]
- **Processed results:** Generated by running these scripts on the raw data

## Important Disclaimer

**This code is designed exclusively for our specific dataset and experimental paradigm.** 

- ❌ **Not a general-purpose analysis toolbox**
- ❌ **Not designed for adaptation to other datasets**  
- ❌ **Hardcoded for our exact data structure and format**
- ✅ **Reproduces our exact published results**
- ✅ **Provides full transparency of our methods**
- ✅ **Enables verification of our findings**

If you're working with different auditory cortex data, you would need to substantially modify these scripts or develop your own analysis approach.

## Citation

If you use this code, please cite:

```bibtex
@article{yaron2025peon,
  title={Auditory cortex neurons that encode negative prediction errors respond to omissions of sounds in a predictable sequence},
  author={Yaron, Amit and Shiramatsu-Isoguchi, Tomoyo and Kern, Franziska B and Ohki, Kenichi and Takahashi, Hirokazu and Chao, Zenas C},
  journal={PLOS Biology},
  year={2025},
  note={Manuscript ID: PBIOLOGY-D-25-00849R2}
}
```

## License

This project is licensed under the MIT License - see [LICENSE](LICENSE) file for details.

## Support & Contact

**For reproducing our paper results:**
- **Issues with our code:** Use GitHub Issues  
- **Questions about our methods:** Contact Amit Yaron at amityar@gmail.com
- **Data access problems:** Check Zenodo DOI above

**Please note:** We provide support for reproducing our published results, but cannot assist with adapting this code for other datasets or experimental paradigms.

## Contributors

- **Amit Yaron** - Primary developer and analysis design
- **Zenas C. Chao** - Senior supervision and methodology

---

**Keywords:** PEON analysis, omission neurons, predictive coding, auditory cortex, reproducible research, Yaron et al 2025
