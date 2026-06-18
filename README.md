## Modeling Psychological Time Series with Multilevel Hidden Markov Models: A Numerical Evaluation

This repository contains the code for a numerical evaluation of the multilevel Hidden Markov Model
with the R-package [mHMMbayes](https://cran.r-project.org/web/packages/mHMMbayes/index.html). 

This repository relates to the paper:

[Aarts, E., & Haslbeck, J. M. B. (2025, July 16). Modeling Psychological Time Series with Multilevel Hidden Markov Models: A Numerical Evaluation. Doi: 10.31234/osf.io/b5mxk_v2](https://osf.io/preprints/psyarxiv/b5mxk_v2)

This repository contains the following files and folders:

- `1_Sim_files_cluster` contains the files to generate 100 independent intensive longitudinal 
datasets from the multilevel HMM for each simulation scenario with given number of states M, pairwise DKL, 
and number of variables, persons, and time points, using the R-package `mHMMbayes`. Within the same files, we fit a 1-, 2-, 3-, 
and 4-state multilevel HMM with the R-package `mHMMbayes`. Files with `p2_` and `p4_`in the name pertain to code for scenarios that were run in 2 or 4 parts in the cluster to save computation time. 
Files that contain `_conv` in the file name pertain to code to run extra chains to inspect convergence for a subset of the 100 simulation iterations. 
For varying *true* number of states, separate generating files are used. Varying the pairwise DKL, 
and number of variables, persons, and time points is done by reading line by line from the file `sim_scenarios_V4.csv`, 
where the columns contain the number of persons $n$, the number of time points $n_t$,the number of variables $ndep$, 
and the pairwise DKL $KLD$. 
- `2_Post-processing` The files `Post-processing j style *st.R` (where \* takes on `1st`, `2st`, or `3st`) contain the code to process the raw output simulation result files from the cluster into `extracted_results_*st.RDS` (step 1; where \* takes on `1st`, `2st`, or `3st`), and further convert these into files containing the relevant performance metrics only, which are stored in the folder `3_results_tables` 
- `3_results_tables` Contains files storing the relevant performance metrics for each simulation scenario and simulation run. See the readme file in this folder for details. 
- `4_plotting results` Contains the code used to visualize the simulation results using the processed output in the folder `3_results_tables`.
- `5_explore_KLD` Contains the code used to calculate the pairwise KLDs for the empirical datasets and visualize (contained in supplement of the paper)
- `6_Sensitivity` Contains the the files related to the sensitivity analysis to the specification of hyper-prior parameters. This includes the 
files used to generate 100 independent intensive longitudinal datasets, post processing files from the raw output simulation result files from the cluster, and
files storing the relevant performance metrics for each simulation scenario and simulation run.

Note that due to the very large size of the output files of the simulation study, the raw output files, and partly processed raw output files (`extracted_results_*st.RDS`, `extracted_convergence_*st.RDS`, and `Convergence_GMStat.RDS`; where \* takes on `1st`, `2st`, or `3st`) are not included in this repository. 
