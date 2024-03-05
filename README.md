# GEn Analysis
SBS GEn analysis framework. This was initially copied over from Provakar Datta's framework for GMn but has since diverged significantly.

# Configuration files

The configuration files are located in the **config** directory. These tells the analysis scripts all the relevant information about a specific run, i.e. beam energy, kinematic angles, raw data location.

# Analyze Raw Data

The main analysis scripts are located in the directory **scripts/replay**. The main script for production data is:
```shell
QuasiElastic_ana.C
```
The main script for simulation data is:
```shell
QuasiElastic_sim_ana.C
```

These files will output a condensed analyzed root file that is much less data than the raw root file. This scripts can take hours to run if you are analyzing hundreds of millions of events

# Asymmetry Analysis

The asymmetry analysis scripts are located in the directory **scripts/analysis**. These scripts are for GEN physics analysis and focus on making elastic cuts and obtaining asymmetry yields. There are many scripts and all of them have header comments explaining what they do. The main script calculating all the asymmetries is
```shell
Asymmetry_yields.C
```