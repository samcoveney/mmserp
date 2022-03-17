# README

You will need openCARP installed. As for Python, dependancies will reveal themselves when attempting to run code.


## Installation

Firstly, compile the stan programs in the directory "stanmodels" by calling `python <filename>.py`.
Then move the resulting <filename>.pkl files to the mmserp/data directory.

Then install as with `python setup.py install`


## Mesh to data storage

To keep everything organized, HDF5 files are used.

The first task is to transfer a simulation ready mesh to a new HDF5 file.

Take a look at the script `mmserp_meshToHDF5`. It assumes that the mesh is stored in some specific files, that can be inferred by reading the code in the script. When you have your mesh stored in this way, call this script.

TODO: remove LGE dependancy in the storage, not needed

TODO: how on earth to extract the data from these? https://zenodo.org/record/5801337#


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


