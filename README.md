# Non-Gaussian Simultaneous Autoregressive Models with Missing Data


This repository contains R code for the simulation studies and real-world examples presented in the manuscript:

- **Simulation1: SEM under MAR n=625, with **75%** missing data (Sim_MAR)**
- **Simulation2: SEM under MNAR n=625, with **75%** missing data (Sim_MNAR)**
- **Real example1:SEM under MAR with **75%** missing data (Real_MAR)**
- **Real example2:SEM under MNAR with **75%** missing data (Real_MNAR)**

## Folder Structure

Each example is stored in a separate folder, which contains the following sub-folders:

- **source.R**: Contains the R code to run the VB algorithms.
- **implement.R**: Contains the R code to implement all algorithms, and generate plots.
- **SEM_MNAR_HMC.stan**: Contains the stan code to run the HMC algorithms (Only for the simulations).

## Running the Scripts

To reproduce results in the manuscript, for example, the **Simulation1: SEM under MAR n=625, with **75%** missing data**:

1. **Download the `Sim_MAR` folder.**
2. **Set the working directory** to this folder.
3. **Run** the `implement.R` script.

Similarly, for other examples, set the working directory to the respective folder before running the corresponding `implement.R` script.

<!--### **Reproducing Pre-Saved Results**

To generate plots and output using pre-saved data:

- Set `rerun_vb` and `rerun_hmc` in the `*_main.R` script to `FALSE`. The script will load results automatically.
- To re-run the VB and HMC algorithms from scratch, set `rerun_vb` and `rerun_hmc` to `TRUE`.

### **Supplementary Sections**

- **Section S3 (Variance Testing)**: Run `var_test_*.R` files in `Logistic/var_test` and `Polypharmacy/var_test` folders.
- **Section S4 (Repeated Simulations)**: Run `*_multi_sims.R` files in `1_Linear/multi_sims/`, `2_Logistic/multi_sims/`, and `5_Poisson/multi_sims/` folders.

The flag `use_tempering` (default: `TRUE`) enables the damped version of VB, as used in the paper.-->

## RStudio and Package Requirements

### **R Version Compatibility**

- The code was tested on **R version 4.1** and **RStan version 2.21**.
- It is also compatible with **R version 4.3** and **RStan version 2.26**.
- Ensure your R installation is configured to compile C++ before installing RStan.

### **Required R Packages**

Install dependencies using:

```r
install.packages(c( "rstan","Matrix","coda","mvnfast","patchwork","vctrs","tidyr","igraph", "ggplot2", "MASS", "spdep","tictoc" ,"mvtnorm", "dplyr","reshape2","spatialreg"))
```


For detailed installation instructions and system requirements, refer to the respective package documentation.
