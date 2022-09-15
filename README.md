# README

[![DOI](https://zenodo.org/badge/445878250.svg)](https://zenodo.org/badge/latestdoi/445878250)


## Installation

You will need openCARP installed.

Firstly, compile the stan programs in the directory "stanmodels" by calling `python <filename>.py`.
Then move the resulting <filename>.pkl files to the mmserp/data directory.

Then install with `python setup.py install` from the top-level folder.

Dependancies should reveal themselves when attempting to run code.


## Using this code

For convenience, this project is organized via scripts, which should be accesible from the Python environment after installation.
Most scripts can be called as `mmserp_scriptName --help` to see the available input arguments.
The rest of this README explains how to use these scripts.


## Mesh to data storage

HDF5 files are use to keep everything organized.

* `mmserp_meshToHDF5` - save a mesh into a new HDF5 file
* `mmserp_browseHDF5` - graphically browse the HDF5 file 
* `mmserp_duplicateHDF5` - create a new HDF5 file, containing only mesh related data, from an existing HDF5 file (requires that the eigenproblem has been solved)


## Solving the Laplacian eigenproblem

The calibration utilizes [Gaussian Process Manifold Interpolation (GPMI)](https://royalsocietypublishing.org/doi/10.1098/rsta.2019.0345)

* `mmserp_viewMesh` - plot the mesh colored by Universal Atrial Coordinates (UACs)
* `mmserp_decimateMesh` - create a lower resolution mesh, and solve Laplacian eigenproblem 
    * use `--loc`, `--num`, and `--runs` to determine vertices of lower resolution mesh
    * use `--tri` to triangulate these vertices
    * use `--holes`, `--layers`, `--eigs`, `--numeigs` to solve eigenproblem
* `mmserp_viewEigs` - plot the mesh colored by a specified eigenvector


## Running atrial simulations

* `mmserp_generateFields` - generate and store random spatially-correlated parameter fields
* `mmserp_createStimulus` - define a set of vertices to use for stimulus in simulations, and save them with a specific name
* `mmserp_createCARPfiles` - create a directory of files for running the CARP simulations
* `runCARP.sh` - run simulations from the directory created above
* `mmserp_getSimResults` - get simulation results from the directory created above


## Calibrating parameter fields

* `mmserp_inference` - perform inference using a virtual experimental design
    * use `--params` to specify the name of parameters as specified in the HDF5 file
    * use `--obs` and `--num` to create and save a design of observation locations
    * use `--deltaS2` to define the resolution of the S1S2(S3) pacing
    * use `--inference` and `--basis` to do MCMC inference of the parameter fields
    * use `--plotGroundTruth`, `--plotAPD` (with `--pacesite`), and `--plotInferred` for plotting


