# Defining a Biomass Limit Reference Point Under Time-Varying Productivity, A Case Study of Atlantic Cod in the Southern Gulf of St. Lawrence, Canada  
*François Turcotte, Jenni L. McDermid, Daniel Ricard, Douglas P. Swain*  

---

## Abstract
Many commonly used fisheries reference points assume stable population processes. However, time-varying productivity is increasingly observed among fish stocks. Atlantic cod (*Gadus morhua*) in the southern Gulf of St. Lawrence (sGSL) collapsed in the 1990s and has not recovered, with non-stationarity in recruitment, growth, maturity, and natural mortality preventing equilibrium conditions. In Canada, the biomass limit reference point (LRP) defines the threshold below which a stock is considered to experience serious harm to its productivity and requires rebuilding. This repository provides data and code used to evaluate candidate LRPs for sGSL cod while accounting for non-stationarity in productivity.  

---

## Repository Contents  
- **`SCA_DATA.RDS`**  
  Annual data including:  
  - Year  
  - Weight-at-age (`waa`)  
  - Maturity-at-age (`mat`)  
  - Total fisheries catch (`catch`)  
  - Statistical catch-at-age model estimates:  
    - Fishery selectivity (`sel`)  
    - Natural mortality (`M`)  
    - Spawning stock biomass (`ssb`)  
    - Number of age-2 fish (`recruits`)  
    - Biomass at ages 2–11 (`b2to11`)  

- **`allee.R`**  
  R script to:  
  - Read data  
  - Calculate production  
  - Perform segmented and linear regressions of production rate vs. biomass  

- **`stock_rec_b0_msy.R`**  
  R script to:  
  - Read data  
  - Estimate stock–recruit relationships  
  - Derive static and dynamic candidate biomass LRPs  

> **Note:** Scripts can be run independently and in any order.  

---

## How to Use  
1. Clone or download this repository  
2. Open R and load the required packages (see script headers for dependencies)  
3. Run the scripts to reproduce analyses:  
   - `allee.R` for production–biomass relationships  
   - `stock_rec_b0_msy.R` for stock–recruit models and LRP estimation  

---

## Citation  
If you use this repository, please cite:  
**Turcotte, F., McDermid, J.L., Ricard, D., & Swain, D.P. (2025). Defining a biomass limit reference point under time-varying productivity: a case study of Atlantic cod in the Southern Gulf of St. Lawrence, Canada. *Canadian Journal of Fisheries and Aquatic Sciences*.**  
