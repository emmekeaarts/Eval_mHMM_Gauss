# jonashaslbeck@protonmail.com; June 19th, 2026

# --------------------------------------------------------
# ---------- What is happening here? ---------------------
# --------------------------------------------------------

# Taking Emmeke's output of the prior sensitivity analysis 
# and making plots in the same style as in the rest of the paper


# Collecting notes for Emmeke:
# 1) The postprocessing file is not reproducible; all the output files are not in the "data" folder, which does not exist
# 2) The Readme.rmd does not fully correspond to the R-code/files, for example it mentions a "emission_performance_Sensitivity.RDS" which does not exist
# 3) 

# --------------------------------------------------------
# ---------- Load Packages -------------------------------
# --------------------------------------------------------


# --------------------------------------------------------
# ---------- Load Preprocessed Output --------------------
# --------------------------------------------------------

# Label Switching
Label_switch_proxy_Sensitivity <- readRDS("6_Sensitivity/result_tables/Label_switch_proxy_Sensitivity.RDS")

# Emission Means
Performance_emission_Sensitivity <- readRDS(file = "6_Sensitivity/result_tables/Performance_emission_Sensitivity.RDS")
extracted_results_Sensitivity <- readRDS("Extracted_results/extracted_results_Sensitivity.RDS")

# Transition Probabilities


# State Decoding



# --------------------------------------------------------
# ---------- More Code of Emmeke -------------------------
# --------------------------------------------------------



# --------------------------------------------------------
# ---------- Plotting: Label Switching -------------------
# --------------------------------------------------------


# --------------------------------------------------------
# ---------- Plotting: Emission Means --------------------
# --------------------------------------------------------


# --------------------------------------------------------
# ---------- Plotting: Transition Probabilities ----------
# --------------------------------------------------------



# --------------------------------------------------------
# ---------- Plotting: State Decoding --------------------
# --------------------------------------------------------









