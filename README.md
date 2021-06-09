# HCVspread

## A multi-scale mathematical model to analyze HCV spread kinetics

The provided code defines a multi-scale mathematical model to analyze viral spread within two-dimensional in vitro cell culture systems following the spatio-temporal dynamics of infection. It relates to the following publications with details on the model structure and implementation provided within the articles:
* Kumberger P, Durso-Cain K, Uprichard SL, Dahari H, Graw F: **“Accounting for Space-Quantification of Cell-To-Cell Transmission Kinetics Using Virus Dynamics Models”**. *Viruses* 2018;10(4). doi: 10.3390/v10040200.
* Durso-Cain K, Kumberger P, Schälte Y, Fink T, Dahari H, Hasenauer J, Uprichard SL, Graw F: **”HCV spread kinetics reveal varying contributions of transmission modes to infection dynamics”**.

## Description:
The simulation environment considers the proliferation and infection of individual cells, as well as viral replication, production and diffusion of viral particles within a 2D monolayer of cells. It comprises a hybrid stochastic-deterministic agent-based model with intracellular viral replication and export kinetics following dynamics described by ordinary differential equations, while intercellular behavior, i.e. viral transmission and cell infection are defined by stochastic events. Cells are modeled as hexagonal shapes and distributed across a predefined grid. The code contains the basic structure of the model, neglecting additional components, such as cell death and antiviral effects that are easily added in other versions of the model. The model has been applied to analyze viral spread dynamics within HCV spread culture systems, as well as to general aspects of the spatio-temporal dynamics of viral spread among cells. All details of the model and its implementation are provided in the aforementioned publications above (Kumberger *et al.* Viruses 2018, Durso-Cain *et al.* 2021).

The model has been developed and implemented by Peter Kumberger and Frederik Graw. It has been adapted for high-performance computing approaches and combined with the parameter inference framework pyABC (Klinger E, Rickert D, Hasenauer J. pyABC: distributed, likelihood-free inference. *Bioinformatics* 2018;34(20):3591-3.doi:10.1093/bioinformatics/bty361; https://pyabc.readthedocs.io/en/ ) in order to apply the model to experimental data on HCV foci spread assays (see also Durso-Cain *et al.* 2021).

## Implementation & Requirements:
The simulation environment is written in C++ and R Version 3.3.3, with no specific requirements with regard to R packages. There is no specific installation needed to run the simulations. Simulations are started in R with several processes running in C++ in the background. The following files build the structure of the code: 

* **rand_generator.h:** File with functions for random number generation as used in the C++ files
* **HCVspread_parameter.h:** File containing the parameters used to run the simulations. Will be automatically generated when running the simulations (*HCVspread_rscript.R*) out of the table *HCVspread_overview.csv*.
* **HCVspread_corefile.h:** File containing the core-functions of the simulation environment defining cell types and viral grid, as well as the main processes (e.g. cell proliferation, viral replication, cell infection by cell-free and cell-to-cell transmission, etc.). Calls the files *HCVspread_parameter.h* and *rand_generator.h*.
* **HCVspread_mainfile.cc:** File containing the actual function to simulate viral spread as it is called from the R function. The function calls the file *HCVspread_corefile.h*.
* **HCVspread_rscript.R:** Actual R-file calling the simulation function and running the simulation.
* **HCVspread_overview.csv:** File containing the individual parameters for each simulation run. The file is called by *HCVspread_rscript.R* to generate the parameter file *HCVspread_parameter.h*.

Some of the paths within the scripts might have to be adapted to the specific user requirements. Output is provided as .csv-files with the different tables containing the following information:

* **HEPstat_xx.csv:** Table with the numbers of each cell type at each time point specified for saving.
* **HEPinf_xx.csv:** Table with each line containing the Cell IDs of all infected cells at each time point specified for saving
* **HEPinfectious_xx.csv:** Table with each line containing the Cell IDs of all infectious cells at each time point specified for saving
* **HEPInfectClusterType_xx.csv:**
* **Empty_xx.csv:**
* **VIRALGRID_xx.csv:**
