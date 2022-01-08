# README

This directory contains a number of scripts for taking a UAC mesh, creating simulaiton files, and extracting simulation results.


## Processing UAC mesh data

* UAC_to_hdf5.py - run this script from a UAC mesh directory to create a HDF5 file for this mesh (saved in group "sim_mesh").
* view_UAC_mesh.py - run this script from anywhere, supplying the HDF5 filename, to view UAC and LGE. 
* decimate.py - run this script from anywhere, supplying the HDF5 filename, to create an interpolation resolution mesh (saved in group "mesh").
* duplicate.py - copy the mesh data from a HDF5 file into a new HDF5 file.


## Creating simulation files

* create_stimulus_file.py - generates a stimulus file (saved in /tmp) for a specific mesh at UACs specified by user.
* upsample_params.py - upsamples the predicted parameters to simulation resolution (save in /tmp/sim_mms_params.npy).
* create_carp_files.py - create a directory of files for running the CARP simulations.


## Simulation 

* runCARP.sh - run the CARP simulations
* get_sim_results.py - retrive simulation results (LAT, APD, |CV|) at interpolation vertex resolution, and save in the HDF5.


## Miscellaneous
* turbocmap.py - contains the turbo colormap (better version of jet), used in a few of these scripts.


