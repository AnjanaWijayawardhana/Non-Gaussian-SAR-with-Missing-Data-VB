# Non-Gaussian Simultaneous Autoregressive Models with Missing Data


This repository contains R code for the simulation studies and real-world examples presented in the manuscript:

---

##  Folder Structure

Simulation studies and the real-world example are organized into separate folders. Each contains:


- **`Source.R`** – Contains the core R code for running the algorithms.  
- **`Implement.R`** – Implements and runs all algorithms for the selected scenario.

---

##  Running the Scripts

### Simulation Studies

To reproduce results from the simulation section of the manuscript (e.g., *Fit SEM under MAR with n = 5,041 and 90% missing data for N=250 datasets*):

1. **Download** the `Simulations` folder.
2. **Set the working directory** to this folder in R or RStudio.
3. **Open** the `Implement.R` script.
4. **Specify**:
   - the model type  
   - the number of observations (`n`)  
   - the number of simulations (`N`)  
   - the percentage of missing data  
5. **Run** the script to perform the simulation.

> To reproduce other simulation results, repeat the steps using the appropriate settings.


###  Real-World Example

To reproduce results from the real-world application (e.g., *Fit SEM to the Lucas County house price data under MAR with 90% missing data*):

1. **Download** the `Real applications` folder.
2. **Set the working directory** to this folder in R or RStudio.
3. **Open** the `Implement.R` script.
4. **Specify** the model type and the missing data percentage.
5. **Run** the script to perform the analysis.


---

## Environment Requirements

###  R Version

- Tested on **R 4.3.1**  
- Also compatible with **R 4.4.2**

### Required Packages

Install the necessary packages by running:

```r
install.packages(c("numDeriv", "Matrix", "spatialreg", "spData", "spdep", "tictoc", "igraph"))
tailed installation instructions and system requirements, refer to the respective package documentation.
