This repository features all the scripts to reproduce the figures and run the models that were produced for:
van Engelen et al. (2019), A three-dimensional palaeohydrogeological reconstruction of the groundwater salinity distribution in the Nile Delta Aquifer
DOI: 10.5194/hess-2019-151

To run the models, use the scripts in the folder "process".
The jobscript Run_Model can be called to start the models, but you will likely need to update the paths in this script.
The following software should be installed on the cluster (running on Linux):
-SLURM as a workload manager
-iMOD-SEAWAT (or iMOD-WQ)
-GNU Parallel
-Miniconda (create python 3.6 environment with the provided .yaml file)
-7zip (optional for backups of modelinput)

In the folder path_management there are .csv files where paths to the folders with modeloutput should be provided.

The scripts analyze_* are used to convert the big netcdf files to smaller files that can be used to plot, either as 1D .nc files or .csv files.
The scripts plot_* are for plotting these smaller datafiles.

SeabornFig2Grid is a utility class that is used as a hack to include seaborn PairGrids in gridspecs.

The version included in the folder might work, but if not:
To get iMOD-SEAWAT executables or to compile iMOD-SEAWAT yourself from the source code, 
either wait for the official iMOD-WQ v5 release (planned in November 2019) and follow the instructions here:
https://oss.deltares.nl/web/imod/deltares-executables-of-imod
or contact imod.support@deltares.nl and ask for a beta version of iMOD-WQ v5