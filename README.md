# Bayesian calibration of ice sheet models
This repository contains code to perform Bayesian calibration of an ice sheet model ensemble. Currently, the code calibrates one particular ensemble of the Greenland Ice Sheet, created using the Ice-sheet and Sea-level System Model (ISSM). Here is a rough outline of the steps involved:
1. Preprocess input data (see the directory GrIS_committedSLR_preprocess/)
1. [Optional] Convert netcdf files to zarr, if running on AWS cloud
1. Perform the calibration (see the directory GrIS_committedSLR_calibration/)
