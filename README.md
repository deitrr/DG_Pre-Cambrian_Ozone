[![DOI](https://zenodo.org/badge/547565921.svg)](https://zenodo.org/badge/latestdoi/547565921)


# DG_Pre-Cambrian_Ozone
Scripts for plots and analysis of Deitrick &amp; Goldblatt pre-Cambrian ozone paper [https://doi.org/10.5194/cp-19-1201-2023](https://doi.org/10.5194/cp-19-1201-2023)

The simulation output necessary to plot all figures is located here: [(https://doi.org/10.20383/103.0648)](https://doi.org/10.20383/103.0648)

## Requirements

This code requires Python 3 as well as the following Python libraries: numpy, scipy, matplotlib, netCDF4, climt, cdo, basemap

### Note about installing climt in Mac OS

Depending on the sub-version of Python 3, installing via `pip install climt` often fails. I've found that I need to install and set the definitions of GCC and GFortran explicitly, like so:
```
$ brew install gcc
$ export CC=gcc-x
$ export FC=gfortran-x
```
where you replace the "x" with whatever version number of gcc brew installs. Then `pip install climt` works for me. (The issue seems to be that Mac OS interprets "GCC" as Apple Clang unless told otherwise)

## Post-processing the simulation output

First, make sure you have downloaded the repository of simulation output (above). Make note of the path where this is located. The directory "DG_Pre-Cambrian_O3_archive" contains two sub-directories, "CAM_input" and "CAM_output". The former contains files used to define the ozone in CAM and the initial condition file, which is only included here for plotting purposes. "CAM_output" contains all simulation output related to the atmosphere for years 31-60, in merged netCDF4 files. These are organized in sub-directories with each simulation name. This output data needs to be averaged and post-processed in a number of ways to be useful for plotting.

Navigate to the directory "CAM_post_process". The post-processing is done using the script `post_process.py`:
```
$ python post_process.py <path_to_archive>/DG_Pre-Cambrian_O3_archive/CAM_output
```
This will process all simulations in the "CAM_output" directory (make sure there are no other files in this directory!). This will take a few hours to run. If the process gets interrupted, you can re-run post_process.py the same way as above and it will detect the existing post-processed files and skip them. You can also run the script for only one simulation by adding the option `-s <simulation_name>`. Use the `-w` option to overwrite all the post-processed outputs, if needed.

`post_process.py` will output post-processed files in sub-directories with the simulation name within directory you execute it from ("CAM_post_process", most likely). The total size of all of these files is about 500 MB. 

## Plotting 

Each figure script in the main directory can be run with
```
$ python figure_x.py
```
where "x" is the figure number. Some of them (`figure_1.py`, `figure_3.py`, and `figure_5_to_7.py`) require data from "DG_Pre-Cambrian_O3_archive/CAM_input/", so you will need to edit those scripts with the correct path to the simulation archive. 

## Credits

Much of the post-processing is done using the [Python implementation of CDO](https://code.mpimet.mpg.de/projects/cdo/wiki/Cdo%7Brbpy%7D).

The Eulerian mean-stream function is calculated using Equation 1 of [Pauluis et al. 2008](https://doi.org/10.1126/science.1159649).

The LTS and EIS are calculated using definitions in [Wood & Bretherton 2006](https://doi.org/10.1175/JCLI3988.1). The EIS is dependent on the height of the LCL. For this we use Equation 24 of [Lawrence 2005](https://doi.org/10.1175/BAMS-86-2-225).

The Python package `climt` is described in [Monteiro et al. 2018](https://doi.org/10.5194/gmd-11-3781-2018). Documentation is available [here](https://github.com/CliMT/climt) and [here](https://climt.readthedocs.io/en/latest/). Documentation specifically for RRTMG is [here](http://rtweb.aer.com/rrtm_frame.html).
