# DG_Pre-Cambrian_Ozone
Scripts for plots and analysis of Deitrick &amp; Goldblatt pre-Cambrian ozone paper (I'll update this paper info when the paper is published)

The simulation output necessary to plot all figures is located here: (I'll update with link when available)

## Requirements

This code requires Python 3 as well as the following Python libraries: numpy, scipy, matplotlib, netCDF4, climt, cdo, basemap

### Note about installing climt in Mac OS

Depending on the subversion of Python 3, installing via `pip install climt` often fails. I've found that I need to install and set the definitions of GCC and GFortran explicitly, like so:
```
$ brew install gcc
$ export CC=gcc-x
$ export FC=gfortran-x
```
where you replace the "x" with whatever version number of gcc brew installs. Then `pip install climt` works for me. (The issue seems to be that Mac OS interprets "GCC" as Apple Clang unless told otherwise)
