# Non-Gaussian Simultaneous Autoregressive Models with Missing Data

This repository contains **R code** for the simulation studies and real-world examples presented in the manuscript.

---

## 📂 Folder Structure

The repository is organized into simulation studies and real-world applications:

- **`Simulation-VB vs HMC`**  
  Code to compare the accuracy of Variational Bayes (VB) algorithms with Hamiltonian Monte Carlo (HMC).

- **`Simulation-compare SEMs`**  
  Code to evaluate the robustness of different Spatial Error Models (SEMs).

- **`Real applications`**  
  Code for fitting different SEMs to the **Lucas-1998-HP dataset**.

Each of these folders contains two subfolders:
- **Full data** – Analysis with complete datasets.  
- **Missing data** – Analysis with incomplete datasets.  

Each subfolder includes:  
- **`Source.R`** – Core R code for running the algorithms.  
- **`Implement.R`** – Implementation code to run the algorithms and generate plots.  
- **`YJ_SEM_Gau_full.stan`** and **`YJ_SEM_Gau_miss.stan`** – Stan models for HMC algorithms (only in the `Simulation-VB vs HMC` subfolders).  

---

## ▶️ Running the Scripts

To reproduce results from the manuscript, for example **comparing VB and HMC methods with full data**:

1. Download the **`Simulation-VB vs HMC`** folder.  
2. Open the **`Full data`** subfolder.  
3. Set your **working directory** to this subfolder in R.  
4. Run the **`Implement.R`** script.  

For other examples, follow the same procedure:  
- Navigate to the relevant folder (full or missing data).  
- Set it as your working directory.  
- Run the `Implement.R` script.  

---

## ⚙️ Environment Requirements

### R Version
- Tested on **R 4.3.1**  
- Compatible with **R 4.4.2**

### Required Packages
Install required packages with:

```r
install.packages(c("numDeriv", "Matrix", "spatialreg", 
                   "spData", "spdep", "tictoc", "igraph"))
