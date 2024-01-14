# MSIM
MSIM is the **MAVIS** fork of the open-source [HARMONI-ELT/HSIM](https://github.com/HARMONI-ELT/HSIM) tool. MSIM has been adapted for simulating observations with MAVIS on the VLT U4. MSIM takes high spectral and spatial resolution input data cubes, encoding physical descriptions of astrophysical sources, and generates mock observed data cubes. The simulations incorporate detailed models of the sky, telescope, instrument, and detectors to produce realistic mock data. If you use MSIM, please consider citing the HSIM paper, [Zieleniewski et al. 2015b](https://doi.org/10.1093/mnras/stv1860).

The purpose of MSIM is to perform in-depth simulations for several key science cases, as part of the design and development of the MAVIS integral field spectrograph. The goal is to be able to simulate the performance of the instrument folding in detailed parameters including the MAVIS AO performance, atmospheric effects and realistic detector statistics.

For updated information on how to run MSIM, you can try the HSIM3 manual ([hsim3.pdf](https://github.com/HARMONI-ELT/HSIM/blob/master/hsim/manual/hsim3.pdf)), though please note that there is no planned long-term compatibility between HSIM and MSIM, so this manual may diverge from MSIM eventually.

## System Requirements
The pipeline is written in Python v3.6. It has been tested on Linux (Ubuntu 18.04 and 20.04), although it should be possible to run on Windows and Mac OSX.

The required Python modules are:
- astropy 4.1
- numpy 1.19.5
- scipy 1.5.4
- matplotlib 3.3.4
- tqdm 4.65.0

The code has been tested with the indicated package version, although more recent releases of these packages are likely to work as well.

## Tips & Tricks ##
We point out here several useful tips that have emerged from our extensive development process.

1. Datacubes can get very large very quickly and fill up the memory even on large computers. The input datacube can be tailored to suit the required needs. For example if one is only interested in the spatial information of an object (e.g. around an emission line) then the cube can be kept very small in the spectral dimension.

2. The flux density units for input datacubes require the usual flux (e.g., erg/s/cm^2/AA) to be divided by the spaxel scale in arcsec. This then gives e.g. erg/s/cm^2/AA/arcsec^2.


## Contact ##

For questions please contact eric.muller@anu.edu.au

Developers: Eric Muller, Jesse Cranney, Trevor Mendel, and of course, the developers of [HARMONI-ELT/HSIM](https://github.com/HARMONI-ELT/HSIM).
