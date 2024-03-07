'''
Calculates PSF, jitter and telescope background and transmission

CHANGELOG:

Version 0.0.0 (2023-10-27)
--------------------------
- Original HARMONI simulator code
Developers: Miguel Pereira Santaella, Laurence Routledge, Simon Zieleniewsk, Sarah Kendrew

Version 1.0.0 (2024-01-10)
--------------------------
- Added a progress bar to the PSF convolution process
- Removed multi-CPU code
Author: Eric Muller (eric.muller@anu.edu.au)

'''
import sys
import logging

import multiprocessing as mp
import signal
import time

import numpy as np
import scipy.constants as sp
from scipy.signal import fftconvolve

from src.config import *

from src.modules.create_psf import create_psf, define_psf
from src.modules.em_model import *

from tqdm import tqdm # CHANGELOG 09-01-2024: Added progress bar


def process_lambda(params, lamb, image, px, py, pback):
    """
    Process the given image with a specified wavelength.

    Args:
        params (list): List of parameters.
        lamb (float): Wavelength value.
        image (ndarray): Input image.
        px (int): Width of the padding plane.
        py (int): Height of the padding plane.
        pback (float): Back emission value.

    Returns:
        tuple: Tuple containing the updated parameters and the convolved image.
    """
    padding_plane = np.zeros((px, py))
    psf = create_psf(lamb)
    padding_plane[psf.shape[0]:psf.shape[0]+image.shape[0], psf.shape[1]:psf.shape[1]+image.shape[1]] = image - pback
    conv_image = fftconvolve(padding_plane, psf) + pback
    return params, conv_image
    
    
counter = None
result_cube = None
bar_str = None
llambs = None

def save_result(results):
    """
    Save the results of the convolution operation.

    Args:
        results (tuple): A tuple containing the result indices and the convolution image.

    Returns:
        None
    """
    global counter, result_cube, bar_str, llambs
    counter = counter + 1
    sys.stdout.write(bar_str.format(int(100.*counter/llambs), counter, llambs) + "\r")
    sys.stdout.flush()

    (i, x0, x1, y0, y1), conv_image = results
    result_cube[i, :, :] = conv_image[y0:y1, x0:x1]


def sim_telescope(input_parameters, cube, back_emission, transmission, ext_lambs, cube_lamb_mask, debug_plots=False, output_file=""):
    ''' Simulates telescope effects
 
    Inputs:
        input_parameters: input dictionary
            exposure_time: Exposure time [s]
            jitter: Residual telescope jitter [mas]
            air_mass: Air mass of the observation
            zenith_seeing: Atmospheric seeing FWHM [arcsec]
            spaxel_scale: spatial pixel (spaxel) scale [mas]
            telescope_temp: Telescope temperature [K]
            ao_mode: MCAO/NOAO//User defined PSF fits file
            ao_star_hmag: AO star H mag
            ao_star_distance: AO star distance [arcsec]
            user_defined_psf: user defined fits PSF
            n_cpus: no. of CPUs to use
            
    
        cube: Input datacube (RA, DEC, lambda)
        back_emission: Input background emission
        transmission: Input transmission
        ext_lambs: extended lambda array [um]
        cube_lamb_mask: mask array to get the lambs of the cube
        
        debug_plots: Produce debug plots
        output_file: File name for debug plots
        
    Outputs:
        cube: Cube including telescope background and emission and convolved with PSF
        back_emission: back_emission including telescope
        PSF: PSF of the first lambda
    '''
    
    # Get telescope reflectivity
    logging.info("Calculating telescope reflectivity")
    
    telescope_reflectivity = load_transmission_curve(ext_lambs, "VLT_mirror_reflectivity.csv", debug_plots, [output_file, "tel"], "telescope transmission")
        
    back_emission *= telescope_reflectivity
    transmission *= telescope_reflectivity
    
    # Get telescope background
    logging.info("Calculating telescope background")
    telescope_background = get_background_emission(ext_lambs, input_parameters["telescope_temp"], 1. - telescope_reflectivity, input_parameters["exposure_time"], debug_plots, [output_file, "tel"], "telescope emission [photons/m$^2$/$\mu$m/arcsec$^2$]")
    back_emission += telescope_background
    
    # Add telescope emission/transmission to the input cube
    tel_reflectivity_cube = telescope_reflectivity[cube_lamb_mask]
    #tel_reflectivity_cube.shape = (np.sum(cube_lamb_mask), 1, 1)
    cube *= tel_reflectivity_cube[:, np.newaxis, np.newaxis]

    tel_background_cube = telescope_background[cube_lamb_mask]
    #tel_background_cube.shape = (np.sum(cube_lamb_mask), 1, 1)
    cube += tel_background_cube[:, np.newaxis, np.newaxis]
    
    # PSF + Jitter + Instrument PSF
    logging.info("Define PSF")
    jitter = input_parameters["jitter"]
    spax = input_parameters["spaxel_scale"]
    
    FWHM_instrument = (config_data["dynamic_instrument_psf"]**2 + config_data["static_instrument_psf"][spax]**2)**0.5
    sigma_instrument = FWHM_instrument/2.35482
    sigma_combined = (jitter**2 + sigma_instrument**2)**0.5
    
    logging.info("Residual telescope jitter = {0:.2f}x{1:.2f} mas".format(*jitter))
    logging.info("Instrument PSF σ = {:.2f} mas".format(sigma_instrument))
    logging.info("-> Combined σ = {0:.2f}x{1:.2f} mas".format(*sigma_combined))
    
    psf_size = config_data["spaxel_scale"][spax].psfsize
    
    define_psf(input_parameters, sigma_combined, psf_size, config_data["spaxel_scale"][spax].psfscale)
    lambs = ext_lambs[cube_lamb_mask]
    
    # padding with back_emission
    padding_x = cube.shape[2] + 2*psf_size
    padding_y = cube.shape[1] + 2*psf_size
    
    padding_back = back_emission[cube_lamb_mask]
    #padding_back.shape = (len(lambs), 1, 1)
    
    conv_side_x = padding_x + psf_size - 1
    conv_side_y = padding_y + psf_size - 1

    # extract convolved image
    x0 = int((conv_side_x - cube.shape[2])/2.)+1
    x1 = x0 + cube.shape[2]
    y0 = int((conv_side_y - cube.shape[1])/2.)+1
    y1 = y0 + cube.shape[1]
    
    logging.info("Convolve with PSF")
    
    global counter, result_cube, bar_str, llambs
    
    bar_str = "[ {:2d}% {:" + str(len(str(len(lambs)))) + "d}/{:" + str(len(str(len(lambs)))) + "d} ]"
    
    # CHANGELOG 29-12-2023: Removed multiprocessing as it was causing issues with the code 
    # ncpus = input_parameters["n_cpus"]
    
    # logging.info(("Using {:d} CPU" + ("s" if ncpus > 1 else "")).format(ncpus))
    # if ncpus > 1:

    #   def init_worker():
    #       signal.signal(signal.SIGINT, signal.SIG_IGN)

    #   pool = mp.Pool(ncpus, init_worker)
        
    #   counter = 0
    #   result_cube = np.zeros_like(cube)
    #   llambs = len(lambs)
        
    #   for j in range(0, len(lambs), ncpus):
    #       result = []
    #       for i in range(j, min([len(lambs), j+ncpus])):
    #           result.append(pool.apply_async(process_lambda, args=((i, x0, x1, y0, y1), lambs[i], cube[i, :, :], padding_x, padding_y, padding_back[i]), callback=save_result))
            
            
    #       try:
    #           while True:
    #               time.sleep(0.5)
    #               if all([r.ready() for r in result]):
    #                   break

    #       except KeyboardInterrupt:
    #           pool.terminate()
    #           pool.join()

        
    #   pool.close()
    #   pool.join()

        
    #   cube = result_cube

    # else:
    for i in tqdm(range(len(lambs))):
        # CHANGELOG 09-01-2024: Added progress bar, commented out below
        # sys.stdout.write(bar_str.format(int(100.*i/len(lambs)), i, len(lambs)) + "\r")
        # sys.stdout.flush()
        # #break
        _, conv_image = process_lambda(i, lambs[i], cube[i, :, :], padding_x, padding_y, padding_back[i])
    
        cube[i, :, :] = conv_image[y0:y1, x0:x1]
    
    central_lambda = (lambs[0] + lambs[-1])*0.5
    return (cube, back_emission, transmission), create_psf(central_lambda), central_lambda


