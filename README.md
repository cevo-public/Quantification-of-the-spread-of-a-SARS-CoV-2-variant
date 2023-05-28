# Quantification of the spread of SARS-CoV-2 variant B.1.1.7 in Switzerland

Here lives the code and data for *Chen, Chaoran, et al. "Quantification of the spread of SARS-CoV-2 variant B.1.1.7 in Switzerland." Epidemics (2021); doi: [10.1016/j.epidem.2021.100480](https://doi.org/10.1016/j.epidem.2021.100480).*

Zenodo archive of this repository: [10.5281/zenodo.7979278](https://doi.org/10.5281/zenodo.7979278)

[Our website](https://cevo-public.github.io/Quantification-of-the-spread-of-a-SARS-CoV-2-variant/) gives a concise description about the data and puts the figures into context. 

**This repository was updated between mid-January and mid-April 2021. For newer information, please visit [CoV-Spectrum](https://cov-spectrum.ethz.ch/).**


## Program

### Overview

The analysis is performed in two steps by different sub-programs. In the first step, taking the overall number of reported cases and the results from sequencing, the `code/analysis_and_plots.py` script estimates the advantage of the variant and creates plots of the estimated proportion and cases through time under different assumptions. The generated plots are stored in `figures/`. Further, the script estimates the number of absolute cases of the variant and writes them into `code/Re/data/`. Then, in the second step, the `code/Re/estimate_Re.R` script takes the computed values from `code/Re/data/` and estimates the Re. The resulting values and plots are stored in `code/Re/figures/`.

### Usage

#### Install Dependencies

- Python and the packages defined in `code/requirements.txt`
- R and the packages `tidyverse, lubridate, viridis, EpiEstim`
- Further, this program requires the [shiny-dailyRe project](https://github.com/covid-19-Re/shiny-dailyRe). Please clone that repository and install its dependencies.


#### Run

1. Run `analysis_and_plots.py`. It does not need further configurations.
2. Open `code/Re/estimate_Re.R` and set the `app_dir` variable to the path of the cloned shiny-dailyRe repository.
3. Run `estimate_Re.R`
