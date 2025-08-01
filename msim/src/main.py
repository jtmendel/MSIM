'''
Code for the bulk framework of the HARMONI
simulator. This code should show the main functions/processes
to move from an input datacube (lambda, y, x) to output cubes

CHANGELOG:

Version 0.0.0 (2023-10-27)
--------------------------
- Original HARMONI simulator code
Developers: Miguel Pereira Santaella, Laurence Routledge, Simon Zieleniewsk, Sarah Kendrew

Version 1.0.0 (2024-01-10)
--------------------------
- Removed MCI options from the simulation configuration parameters.
- Manually set the MCI parameter to false.
- Changed from HARMONI AO options to MAVIS AO options.
- Removed code for unused AO options.
- Removed code for unused spaxel scale checks for HARMONI detectors.
- Changed parameter prefixes from HARMONI to MAVIS.
- Changed some plotting code labels from HARMONI to MAVIS.
- Added some logging info for troubleshooting.
- Remover debug plots from being saved to stop code breaking.
Author: Eric Muller (eric.muller@anu.edu.au)

Version 1.0.1 (2024-01-)
- Removed cross-talk application as we have no IR detectors.
- Removed dedicated on-sky observation code as we have no IR detectors.
- Removed detector systematics fits outputs as we have no IR detectors.
#TODO:
'''
import collections
import glob
import re
import datetime
import multiprocess as mp
import os.path
import logging
import warnings

import numpy as np
from scipy.interpolate import interp1d
from scipy.signal import convolve2d
from astropy.io import fits
import astropy.constants as const

from src.config import *
from src.init_cube import init_cube
from src.sim_sky import sim_sky
from src.sim_telescope import sim_telescope
from src.sim_instrument import sim_instrument
from src.sim_detector import sim_detector, mask_saturated_pixels
from src.modules.adr import apply_adr
from src.modules.rebin import *

## Detector systematics
#from src.sim_detector import make_rn_dist, make_dets, add_detectors
#from src.modules.misc_utils import trim_cube


def main(input_parameters):
    '''
    Inputs:
        Dictionary containing:

        'input_cube': Input high resolution datacube (RA, DEC, lambda)
        'output_dir': Output file directory
        'spaxel_scale': spatial pixel (spaxel) scale [mas]
        'grating': Spectral grating
        'n_exposures': NDIT: No. of exposures
        'exposure_time': Exposure time [s]
        
        'zenith_seeing': Atmospheric seeing FWHM [arcsec]
        'air_mass':  Air mass of the observation
        'ao_mode': AO mode
        'user_defined_psf': User defined PSF
        'ao_star_distance': AO star distance [arcsec]
        'ao_star_hmag': AO star H magnitude
        'moon_illumination': Fractional Moon illumination
        
        'extra_jitter': Residual telescope jitter. 1 or 2 x separated numbers [mas]
        'detector_systematics': Boolean - use detector systematics
        'detector_tmp_path':  Directory to save interim detector files
        'telescope_temp': Telescope temperature [K]
        'adr': Boolean - turn ADR on or off
        'mci': Boolean - Use minimum compliant instrument
        'detector': Near-IR detector performance
        
        'noise_seed': 100,
        'config_file': configuration file
        'debug': keep debug plots and all noise outputs
        'n_cpus': Number of processes
        'version': MSIM version
    Outputs:

    '''
    warnings.filterwarnings("ignore", module="astropy")
    warnings.filterwarnings("ignore", module="matplotlib")
    debug_plots = True

    Conf = collections.namedtuple('Conf', 'name, header, value')
    simulation_conf = [
            Conf('MSIM Version', 'MSM_VER', 'version'),
            Conf('Config file', 'MSM_CFG', 'config_file'),
            Conf('Filename', 'MSM_FILE', 'input_cube'),
            Conf('Output dir', 'MSM_OUTD', 'output_dir'),
            Conf('Exposure time', 'MSM_EXP', 'exposure_time'),
            Conf('Number of exposures', 'MSM_NEXP', 'n_exposures'),
            Conf('Spaxel scale', 'MSM_SPAX', 'spaxel_scale'),
            Conf('Grating', 'MSM_GRAT', 'grating'),
            Conf('Zenith seeing', 'MSM_SEEI', 'zenith_seeing'),
            Conf('Air Mass', 'MSM_AIRM', 'air_mass'),
            Conf('Extra jitter', 'MSM_JITT', 'extra_jitter'),
            Conf('Telescope temperature [K]', 'MSM_TEMP', 'telescope_temp'),
            Conf('FPRS temperature [C]', 'MSM_FPRS', 'fprs_temp'),
            Conf('Scattered background [%]', 'MSM_SCSK', 'scattered_sky'),
            Conf('Moon', 'MSM_MOON', 'moon_illumination'),
            Conf('ADR', 'MSM_ADR', 'adr'),
            Conf('Seed', 'MSM_SEED', 'noise_seed'),
            Conf('AO', 'MSM_AO', 'ao_mode'),
            Conf('No. of processes', 'MSM_NPRC', 'n_cpus'), 
            Conf('Internal spectral oversampling factor', 'MSM_SPES', 'spectral_sampling'),
            Conf('Internal spatial oversampling factor', 'MSM_SPAS', 'spatial_sampling'),
            ]

    def str2bool(var):
        if type(input_parameters[var]) is str:
            input_parameters[var] = input_parameters[var].upper() == "TRUE"

    # change type to bool
    str2bool("debug")
    str2bool("adr")
    # str2bool("mci") # CHANGELOG 29-11-2023: Removed the mci option, set to false
    #input_parameters["mci"] = False
    #det_switch = False # CHANGELOG 10-01-2024: Make sure the detector systematics are off
    # str2bool("detector_systematics") # CHANGELOG 10-01-2024: see above line. This option has been removed


    # Init logger
    base_name = os.path.splitext(os.path.basename(input_parameters['input_cube']))[0]
    base_filename = os.path.join(input_parameters['output_dir'], base_name)

    logfile = base_filename + ".log"
    open(logfile, 'w').close()

    # log to a file
    logging.basicConfig(filename=logfile, level=logging.INFO, format='%(asctime)s  %(levelname)s  %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    # and also to the terminal
    logger = logging.getLogger()
    std = logging.StreamHandler()
    std.setFormatter(MSIMFormatter())
    std.setLevel(logging.INFO)
    logger.addHandler(std)

    msimlog = MSIMLoggingHandler()
    logger.addHandler(msimlog)


    # CHANGLOG 29-11-2023: Removed the minimum compliant params
    # if input_parameters["mci"]:
    #     logging.info("mci: forcing air mass = 1.3 = 40deg")
    #     input_parameters["air_mass"] = 1.3
    #     if input_parameters["ao_mode"] == "LTAO":
    #         input_parameters["ao_mode"] = "MCI_LTAO"
    #     elif input_parameters["ao_mode"] == "SCAO":
    #         input_parameters["ao_mode"] = "MCI_SCAO"
    #     else:
    #         logging.error("No valid AO mode for MCI.")
    #         return

    # CHANGELOG 10-01-2024: Removed this code as we use visible detectors, not IR
    # if input_parameters["detector_systematics"] == True:
    #     simulation_conf.append(Conf('Detectors tmp path', 'HSM_DDIR', 'detector_tmp_path'))

    # CHANGELOG 29-12-2023: Removed the LTAO and HCAO options being appended, added MCAO option
    # if input_parameters["ao_mode"] == "LTAO":
    #     simulation_conf.append(Conf('LTAO star H mag', 'HSM_AOMA', 'ao_star_hmag'))
    #     simulation_conf.append(Conf('LTAO star H mag', 'HSM_AODI', 'ao_star_distance'))
    # elif input_parameters["ao_mode"] == "HCAO":
    #     simulation_conf.append(Conf('HC apodizer', 'HSM_HCAP', 'hc_apodizer'))
    #     simulation_conf.append(Conf('HC mask', 'HSM_HCMK', 'hc_fp_mask'))
    if input_parameters["ao_mode"] == "MCAO":
        #TODO: fill this out and append any necessary options
        pass
    elif input_parameters["ao_mode"] == "User":
        simulation_conf.append(Conf('User defined PSF file', 'MSM_UPSF', 'user_defined_psf'))

    # CHANGELOG 29-12-2023: Removed this check
    # # Check that the 60x60 and 120x60 are only used with the V+R grating
    # if input_parameters["spaxel_scale"] in ["60x60", "120x60"] and input_parameters["grating"] != "V+R":
    #     logging.error(input_parameters["spaxel_scale"] + ' is only available for the V+R grating.')
    #     return

    # CHANGELOG 29-12-2023: Removed this option
    # # Check HCAO configuration
    # if input_parameters["ao_mode"] == "HCAO":

    #     logging.warning("HCAO mode is experimental. Please use with caution the results.")

    #     if input_parameters["spaxel_scale"] !=  "4x4":
    #         logging.error("4x4 spaxel scale must be used for the HCAO mode.")
    #         return
    #     if input_parameters["grating"] in ["V+R", "Iz", "z-high"]:
    #         logging.error("V+R, Iz, and z-high gratings are not compatible with the HCAO mode.")
    #         return

    #     if input_parameters["exposure_time"]*input_parameters["n_exposures"] > 10.:
    #         logging.warning("The total exposure time (DIT*NDIT) should be < 10s since field rotation is not simulated by HSIM")

    #     if input_parameters["adr"]:
    #         logging.warning("Disabling standard ADR simulation for HCAO")
    #         input_parameters["adr"] = False

    # Get oversampling factor
    # spectral axis
    if input_parameters["spectral_sampling"] == -1: # Use default oversampling factor
        input_parameters["spectral_sampling"] = config_data["spectral_sampling"]["internal"]
    elif  input_parameters["spectral_sampling"] <= 0:
        logging.error("Spectral sampling must be > 0")
        return
    else:
        if input_parameters["spectral_sampling"] < config_data["spectral_sampling"]["internal"]:
            logging.warning("The selected spectral oversampling (" + str(input_parameters["spectral_sampling"]) + ") is lower than the defaul value: " + str(config_data["spectral_sampling"]["internal"]))

        config_data["spectral_sampling"]["internal"] = input_parameters["spectral_sampling"]

    # Update FPRS temperature
    config_data["MAVIS_FPRS_temp"] = input_parameters["fprs_temp"] # CHANGELOG 10-01-2024: Changed prefix from HARMONI to MAVIS

    # spatial axes
    spax_scale = config_data['spaxel_scale'][input_parameters["spaxel_scale"]]
    if input_parameters["spatial_sampling"] == -1: # Use default oversampling factor
        factor = spax_scale.xscale/spax_scale.psfscale
        input_parameters["spatial_sampling"] = factor
    elif  input_parameters["spatial_sampling"] <= 0:
        logging.error("Spatial sampling must be > 0")
        return
    else:
        psf_fov = spax_scale.psfscale*spax_scale.psfsize
        new_psf_scale = spax_scale.xscale/input_parameters["spatial_sampling"]
        # Update scale info
        spax_scale = SpaxelScaleInfo(spax_scale.xscale, spax_scale.yscale, new_psf_scale, round(psf_fov/new_psf_scale/2.)*2)

        if new_psf_scale > config_data['spaxel_scale'][input_parameters["spaxel_scale"]].psfscale:
            logging.warning("The selected spatial oversampling results in internal " + str(new_psf_scale) + " mas pixels that are larger than the default size of " + str(config_data['spaxel_scale'][input_parameters["spaxel_scale"]].psfscale) + " mas")

        config_data['spaxel_scale'][input_parameters["spaxel_scale"]] = spax_scale

    try:
        res_jitter = input_parameters["extra_jitter"]
        if "x" in res_jitter:
            jitter = np.array([*map(float, res_jitter.split("x"))])
            if len(jitter) != 2:
                logging.error(str(res_jitter) + " is not a valid jitter value.")
                return
        else:
            jitter = np.repeat(float(res_jitter), 2)

        input_parameters["jitter"] = jitter

    except:
        logging.error(str(res_jitter) + " is not a valid jitter value.")
        return


    logging.info("Simulation input parameters:")
    for _ in simulation_conf:
        logging.info(_.value + " = " + str(input_parameters[_.value]))
        if _.value == "config_file":
            logging.info("# start configuration file")

    logging.info("# end configuration file")

    #initialize random number generation
    rng = np.random.default_rng(input_parameters["noise_seed"])
    
    # Read input FITS cube and resample depending on grating and spaxel scale
    # output is in ph/s/m2/um/arcsec2 units
    cube_exp, head, lambs, input_spec_res = init_cube(input_parameters["input_cube"], input_parameters["grating"], input_parameters["spaxel_scale"])

    # Calculate extended lambda range to convolve with the LSF
    LSF_width = int(config_data["spectral_sampling"]["internal"]/2.35482*config_data['LSF_kernel_size'])
    lambs_extended_index = (np.arange(LSF_width)+1)*(lambs[1] - lambs[0])

    lambs_extended = np.concatenate(((lambs[0] - lambs_extended_index)[::-1], lambs, lambs[-1] + lambs_extended_index))
    cube_lamb_mask = np.concatenate(([False]*len(lambs_extended_index), [True]*len(lambs), [False]*len(lambs_extended_index)))

    # calculate the cube in photons for a single exposure
    # in ph/m2/um/arcsec2
    cube_exp *= input_parameters["exposure_time"]
    back_emission = np.zeros(len(lambs_extended)) # Background exposure
    transmission = np.ones(len(lambs_extended)) # Telluric correction

    lambda_data = (lambs_extended, cube_lamb_mask)

    simulated_data = (cube_exp, back_emission, transmission)

    # 1 - Atmosphere: 
    #    - Sky emission (lines + continuum)
    #    - Sky transmission
    #    - Moon
    #    - Atmospheric differential refraction
    simulated_data = sim_sky(input_parameters, *simulated_data, head, *lambda_data, input_spec_res, 
                             debug_plots=debug_plots, output_file=base_filename)

    # 2 - Telescope:
    #    - PSF + Jitter
    #    - Telescope background
    #    - Telescope transmission
    simulated_data, psf_internal, psf_lambda = sim_telescope(input_parameters, *simulated_data, *lambda_data, 
                                                             debug_plots=debug_plots, output_file=base_filename)

    # 3 - Instrument:
    #    - LSF
    #    - Instrument background
    #    - Instrument transmission
    simulated_data, LSF_width_A = sim_instrument(input_parameters, *simulated_data, *lambda_data, input_spec_res, 
                                                 debug_plots=debug_plots, output_file=base_filename)

    cube_exp, back_emission, transmission, fpm_mask = simulated_data

    # 4 - Rebin cube to output spatial and spectral pixel size
    logging.info("Rebin data")
    logging.info("- Output spaxel scale: " + str(input_parameters["spaxel_scale"]))
    z, y, x = cube_exp.shape

    # rebin spatial axes
    scale_x = spax_scale.xscale/spax_scale.psfscale
    scale_y = spax_scale.yscale/spax_scale.psfscale
    out_size_x = int(x/scale_x)
    out_size_y = int(y/scale_y)
    output_cube = np.zeros((z, out_size_y, out_size_x))

    for k in np.arange(0, z):
        output_cube[k, :, :] = frebin2d(cube_exp[k, :, :], (out_size_x, out_size_y))

    if fpm_mask is not None:
        fpm_mask_rebin = frebin2d(fpm_mask, (out_size_x, out_size_y))

    # and update header
    head['CDELT1'] = spax_scale.xscale*np.sign(head['CDELT1'])
    head['CDELT2'] = spax_scale.yscale*np.sign(head['CDELT2'])

    # rebin spectral axis
    scale_z = config_data["spectral_sampling"]["output"]/config_data["spectral_sampling"]["internal"]

    lambs = lambs_extended[cube_lamb_mask]
    new_lamb_per_pix = (lambs[1] - lambs[0])/scale_z # micron

    output_lambs = new_lamb_per_pix*np.arange(int(len(lambs)*scale_z)) + lambs[0]

    logging.info("- Output spectral sampling: {:.2f} AA".format(new_lamb_per_pix*10000.))

    output_cube_spec = rebin_cube_1d(output_lambs, lambs, output_cube)
    output_back_emission = rebin1d(output_lambs, lambs_extended, back_emission)
    output_transmission = rebin1d(output_lambs, lambs_extended, transmission)

    # and update header
    head['CRPIX3'] = 1
    head['CRVAL3'] = output_lambs[0]
    head['CDELT3'] = new_lamb_per_pix
    head['NAXIS3'] = len(output_lambs)


    # correct for ADR
    if input_parameters["adr"] == True:
        logging.info("- Correcting ADR")
        output_cube_spec = apply_adr(output_cube_spec, head, output_lambs, input_parameters["telescope_temp"], input_parameters["air_mass"], correct=True)


    # 5 - Convert to photons per pixel
    spaxel_area = spax_scale.xscale/1000.*spax_scale.yscale/1000. # arcsec2
    channel_width = new_lamb_per_pix # micron

    output_cube_spec = output_cube_spec*spaxel_area*channel_width*config_data["telescope"]["area"]
    output_back_emission = output_back_emission*spaxel_area*channel_width*config_data["telescope"]["area"]
    head['BUNIT'] = "photon"

    # 6 - Detector
    #    - QE
    #    - Dark
    #    - Read noise
    #    - Thermal background
    grating = input_parameters["grating"]

    output_cube_spec, output_back_emission, output_transmission, read_noise, dark_current, thermal_background = \
                sim_detector(input_parameters, output_cube_spec, output_back_emission, output_transmission, output_lambs, \
                    debug_plots=debug_plots, output_file=base_filename)
    head['BUNIT'] = "electron"


    # Generate noiseless outputs
    # - object + background
    output_cube_spec = output_cube_spec.astype(np.float32)

    # - mask saturated pixels
    output_cube_spec, saturated_obj_back = mask_saturated_pixels(output_cube_spec, grating)

    # - background - as a cube the 1D background spectrum
    output_back_emission = output_back_emission.astype(np.float32)
    output_back_emission.shape = (len(output_back_emission), 1, 1)
    output_back_emission_cube = np.zeros_like(output_cube_spec) + output_back_emission

    # - mask saturated pixels
    output_back_emission_cube, saturated_back = mask_saturated_pixels(output_back_emission_cube, grating)

    # - object - back
    output_cube_spec_wo_back = np.subtract(output_cube_spec, output_back_emission_cube)

    # Generate observed outputs
    zero_cube = np.zeros_like(output_cube_spec)

    # A. Object exposure
    DIT = input_parameters["exposure_time"]
    NDIT = input_parameters["n_exposures"]

    # - object exposure with background
    try:
        sim_object_plus_back = rng.poisson(abs(output_cube_spec*NDIT)).astype(np.float32)
    except:
        sim_object_plus_back = rng.normal(loc=abs(output_cube_spec*NDIT), scale=np.sqrt(abs(output_cube_spec*NDIT))).astype(np.float32)

    if np.sum(saturated_obj_back) > 0:
        logging.warning(str(np.sum(saturated_obj_back)) + " pixels are saturated in the obj + back frames")
        sim_object_plus_back[saturated_obj_back] = np.nan

    dark_cube = np.zeros_like(output_cube_spec) + dark_current*NDIT
    thermal_cube = np.zeros_like(output_cube_spec) + thermal_background*NDIT

    # - read noise and dark current for object exposure
    sim_read_noise1 = rng.normal(zero_cube, np.sqrt(NDIT)*read_noise).astype(np.float32)
    sim_dark_current1 = rng.poisson(dark_cube).astype(np.float32)
    sim_thermal1 = rng.poisson(thermal_cube).astype(np.float32)

    # - combine object, read noise and dark current
    sim_total = sim_object_plus_back + sim_read_noise1 + sim_dark_current1 + sim_thermal1

    # B. Background exposure
    #- background only
    try:
        sim_back = rng.poisson(abs(output_back_emission_cube*NDIT)).astype(np.float32)
    except:
        sim_back = rng.normal(loc=abs(output_back_emission_cube*NDIT), scale=np.sqrt(abs(output_back_emission_cube*NDIT))).astype(np.float32)

    if np.sum(saturated_back) > 0:
        logging.warning(str(np.sum(saturated_back)) + " pixels are saturated in the back frames")
        sim_back[saturated_back] = np.nan

    # - read noise and dark current for background exposure
    sim_read_noise2 = rng.normal(zero_cube, np.sqrt(NDIT)*read_noise).astype(np.float32)
    sim_dark_current2 = rng.poisson(dark_cube).astype(np.float32)
    sim_thermal2 = rng.poisson(thermal_cube).astype(np.float32)

    # - combine object, read noise and dark current
    sim_total_only_back = sim_back + sim_read_noise2 + sim_dark_current2 + sim_thermal2

    # C. Calculate reduced cube: object - background exposure
    sim_reduced = sim_total - sim_total_only_back


    logging.info("Pipeline interpolation effects")
    # Convolve the reduced cube with a 1pix FWHM Gaussian to account for the 
    # interpolation during the data reduction
    sigma = 1./2.35482 # pix
    kernel_size = 5

    Gauss2D = lambda x, y: np.exp(-(x**2 + y**2)/(2.*sigma**2))
    xgrid = np.linspace(1, kernel_size, kernel_size) - kernel_size*0.5 - 0.5
    ygrid = np.linspace(1, kernel_size, kernel_size) - kernel_size*0.5 - 0.5
    xx, yy = np.meshgrid(xgrid, ygrid)
    kernel = Gauss2D(xx, yy)
    kernel = kernel/np.sum(kernel)
    z, y, x = sim_reduced.shape
    for pz in range(z):
        sim_reduced[pz, :, :] = convolve2d(sim_reduced[pz, :, :], kernel, mode="same", boundary="fill", fillvalue=0.)

    # Calculate noise cube
    noise_cube_object = abs(output_cube_spec_wo_back*NDIT) # object noise variance
    noise_cube_back = abs(output_back_emission_cube*NDIT) # back noise variance
    noise_cube_read_noise = zero_cube + np.sqrt(NDIT)*read_noise # read noise sigma
    noise_cube_dark = dark_cube # dark noise variance
    noise_cube_thermal = thermal_cube # thermal noise variance

    # CHANGELOG 11-01-2024: Removed the IR detector option
    # if grating != "V+R":
    #     n_observations = 2
    # else:
    n_observations = 1 # no dedicated sky observation

    noise_cube_total = np.sqrt(noise_cube_object + n_observations*noise_cube_back + n_observations*noise_cube_dark 
                               + n_observations*noise_cube_thermal + n_observations*noise_cube_read_noise**2)

    logging.info("Saving output")
    debug_plots = False # CHANGELOG 10-01-2024: Changed debug_plots to false to stop code from breaking
    #TODO: figure out why below code block breaks
    if debug_plots:
        import matplotlib.pyplot as plt

        plt.clf()

        colors = [u'#1f77b4', u'#ff7f0e', u'#2ca02c', u'#d62728', u'#9467bd', u'#8c564b', u'#e377c2', u'#7f7f7f', u'#bcbd22', u'#17becf', u'#2dd42d', u'#eaff00', u'#202020', u'#66f2ff']*2

        fig = plt.figure()
        ax = fig.add_subplot(111)

        total_telescope_sky_em = np.zeros_like(lambs_extended)

        # photon/m2/um/arcsec2 -> W/m2/um/sr
        ph2en_conv_fac = 1./DIT*2.99792458e8/(lambs_extended*1.e-6)*6.62607015e-34*4.25451702962e10

        # a. Sky contribution to the background
        _, tel_tr = np.loadtxt(base_filename + "_tel_tr.txt", unpack=True)

        w, e = np.loadtxt(base_filename + "_sky_em.txt", unpack=True)
        plt.plot(w, tel_tr*e*ph2en_conv_fac, label="sky", color=colors[-1])

        if input_parameters["moon_illumination"] > 0.:
            w, e = np.loadtxt(base_filename + "_moon_em.txt", unpack=True)
            plt.plot(w, tel_tr*e*ph2en_conv_fac, label="Moon", color=colors[6])

        # b. Telescope background emission
        w, e = np.loadtxt(base_filename + "_tel_em.txt", unpack=True)
        plt.plot(w, e*ph2en_conv_fac, label="telescope", color=colors[-2])

        #TODO: change this 
        #if not input_parameters["mci"]:
        # MAVIS parts
        total_instrument_em = np.zeros_like(lambs_extended)
        total_instrument_tr = np.ones_like(lambs_extended)
        mavis_files_em = sorted(glob.glob(base_filename + "_MAVIS_*_em.txt"))

        for mavis_file, color in zip(mavis_files_em, colors):
            # Read part emission
            w, e = np.loadtxt(mavis_file, unpack=True)
            m = re.search('.+MAVIS_(.+)_em.txt', mavis_file)
            # and throughput
            w, t = np.loadtxt(base_filename + "_MAVIS_" + m.group(1) + "_tr.txt", unpack=True)
            plt.plot(w, e/total_instrument_tr*ph2en_conv_fac, label=m.group(1), color=color, ls="--", lw=1.2)
            total_instrument_em = total_instrument_em*t + e
            total_instrument_tr = total_instrument_tr*t

        plt.plot(w, total_instrument_em/total_instrument_tr*ph2en_conv_fac, label="MAVIS total", color="red")
        logging.info("MAVIS emission at input focal plane at {:.4f} um = {:.4e} W/m2/um/sr".format(np.median(w), np.median(total_instrument_em/total_instrument_tr*ph2en_conv_fac)))

        np.savetxt(base_filename + "_total_MAVIS_em.txt", np.c_[w, total_instrument_em], comments="#", header="\n".join([
            'TYPE: Total emission.',
            'Wavelength [um], emission']))

        plt.legend(prop={'size': 6})
        plt.xlabel(r"wavelength [$\mu$m]")
        plt.ylabel(r"back at HRM input [W/m$^2$/$\mu$m/sr]")
        plt.yscale("log")
        plt.savefig(base_filename + "_total_em.pdf")


        ## Transmission

        plt.clf()
        w, e = np.loadtxt(base_filename + "_sky_tr.txt", unpack=True)
        plt.plot(w, e, label="sky", color=colors[-1])
        total_tr_w = w
        total_tr = e

        w, e = np.loadtxt(base_filename + "_tel_tr.txt", unpack=True)
        plt.plot(w, e, label="telescope", color=colors[-2])
        if np.sum(np.abs(total_tr_w - w)) != 0.:
            logging.error('Telescope transmission wavelength error. This should never happen.')
            return
        total_tr *= e

        #TODO: change this
        total_instrument_tr = np.ones_like(total_tr)
        # MAVIS parts
        mavis_files_tr = sorted(glob.glob(base_filename + "_MAVIS_*tr.txt"))
        for mavis_file, color in zip(mavis_files_tr, colors):
            w, e = np.loadtxt(mavis_file, unpack=True)
            m = re.search('.+MAVIS_(.+)_tr.txt', mavis_file)
            plt.plot(w, e, label=m.group(1), color=color, ls="--", lw=1.2)
            total_instrument_tr *= e
            total_tr *= e

        plt.plot(w, total_instrument_tr, label="MAVIS total", color="red") # CHANGELOG 10-01-2024: Changed text from HARMONI to MAVIS
        
        np.savetxt(base_filename + "_total_MAVIS_tr.txt", np.c_[w, total_instrument_tr], comments="#", header="\n".join([ # CHANGELOG 10-01-2024: Changed text from HARMONI to MAVIS
            'TYPE: Total transmission.',
            'Wavelength [um], transmission']))

        # the detector curve has a different wavelength range and spacing
        w, e = np.loadtxt(base_filename + "_det_qe_tr.txt", unpack=True)
        plt.plot(w, e, label="detector", color=colors[5])

        total_trans_interp = interp1d(total_tr_w, total_tr, kind='linear', bounds_error=False, fill_value=0.)
        total_tr_final_w = w
        total_tr_final = total_trans_interp(total_tr_final_w)*e

        plt.plot(total_tr_final_w, total_tr_final, label="Total", color=colors[7])

        plt.legend(prop={'size': 6})
        plt.xlabel(r"wavelength [$\mu$m]")
        plt.ylabel(r"transmission")
        plt.savefig(base_filename + "_total_tr.pdf")

        if input_parameters["debug"] == False:
            list_files = ["sky_tr", "sky_em", "moon_em", "tel_tr", "tel_em", "det_qe_tr"]
            for _ in list_files:
                os.remove(base_filename + "_" + _ + ".txt")
                os.remove(base_filename + "_" + _ + ".pdf")

            for _ in mavis_files_em+mavis_files_tr:
                filename = _.replace(".txt", "")
                os.remove(filename + ".txt")
                os.remove(filename + ".pdf")


    # Noiseless
    outFile_noiseless_object = base_filename + "_noiseless_obj.fits"
    outFile_noiseless_background = base_filename + "_noiseless_back.fits"
    outFile_noiseless_object_plus_back = base_filename + "_noiseless_obj_plus_back.fits"

    # With noise.
    # - Observed Obj + Back + Noise1
    # - Background  Back + Noise2
    # - Reduced  Obj + Back + Noise1 - (Back + Noise2)
    outFile_observed = base_filename + "_observed_obj_plus_back.fits"
    outFile_observed_back = base_filename + "_observed_back.fits"
    outFile_reduced = base_filename + "_reduced.fits"

    # SNR
    # - Obj/Noise
    outFile_SNR = base_filename + "_reduced_SNR.fits"
    outFile_detSNR = base_filename + "_detector_SNR.fits"
    
    # standard deviation cube
    outFile_std = base_filename + "_std.fits"
    outFile_detstd = base_filename + "_det_std.fits"
    
    # detectors
    outFile_alldets = base_filename + "_all_dets.fits"
    outFile_useddets = base_filename + "_used_dets.fits"

    #
    outFile_read_noise = base_filename + "_read_noise.fits"
    outFile_dark = base_filename + "_dark.fits"
    outFile_ddetector_thermal = base_filename + "_det_thermal.fits"

    # Update header
    for _ in simulation_conf:
        head[_.header] = (str(input_parameters[_.value]), _.name)

    head['MSM_TIME'] = str(datetime.datetime.utcnow())

    save_fits_cube(outFile_noiseless_object_plus_back, output_cube_spec*NDIT, "Noiseless O+B", head)
    save_fits_cube(outFile_noiseless_background, output_back_emission_cube*NDIT, "Noiseless B", head)
    save_fits_cube(outFile_noiseless_object, output_cube_spec_wo_back*NDIT, "Noiseless O", head)

    save_fits_cube(outFile_observed, sim_total, "Observed O+B1+Noise1", head)
    save_fits_cube(outFile_observed_back, sim_total_only_back, "Observed B2+Noise2", head)
    save_fits_cube(outFile_reduced, sim_reduced, "Reduced (O+B1+Noise1) - (B2+Noise2)", head)


    save_fits_cube(outFile_SNR, output_cube_spec_wo_back*NDIT/noise_cube_total, "SNR (O-B)/Noise", head)
    save_fits_cube(outFile_std, noise_cube_total, "Noise standard deviation", head)


    # Flux calibration - from electrons to erg/s/cm2/um

    # - Calibration star
    flux_cal_star = 1e-13 # erg/s/cm2/um
    en2ph_conv_fac = (1.98644582e-25 * 1.e7)/(output_lambs*1.e-6*1.e4)
    flux_cal_star_photons = flux_cal_star/en2ph_conv_fac #photon/s/cm2/um

    # Calculate the flux in the PSF core 
    peak_psf = np.max(psf_internal)
    flux_fraction_psf_core = 1. # np.sum(psf_internal[psf_internal > 0.5*peak_psf])

    flux_cal_star_electrons = flux_cal_star_photons*channel_width*config_data["telescope"]["area"]*output_transmission*flux_fraction_psf_core # electron/s 
    factor_calibration = flux_cal_star/flux_cal_star_electrons # erg/s/cm2/um / (electron/s)

    # Reshape factor to match output cube shape
    output_shape = output_cube_spec_wo_back.shape
    factor_calibration = np.repeat(factor_calibration, output_shape[1], axis=0).reshape(output_shape[0], output_shape[1]) 
    factor_calibration = np.repeat(factor_calibration, output_shape[2], axis=0).reshape(output_shape[0], output_shape[1], output_shape[2]) 


    outFile_flux_cal_noiseless = base_filename + "_noiseless_obj_flux_cal.fits"
    outFile_flux_cal_reduced = base_filename + "_reduced_flux_cal.fits"
    outFile_flux_cal_variance = base_filename + "_reduced_variance.fits"


    head['BUNIT'] = "erg/s/cm2/um/arcsec2"
    save_fits_cube(outFile_flux_cal_noiseless, output_cube_spec_wo_back*factor_calibration/(DIT*spaxel_area), "Flux cal Noiseless O", head)
    save_fits_cube(outFile_flux_cal_reduced, sim_reduced*factor_calibration/(NDIT*DIT*spaxel_area), "Flux cal Reduced (O+B1+Noise1) - (B2+Noise2)", head)
    save_fits_cube(outFile_flux_cal_variance, (noise_cube_total*factor_calibration/(NDIT*DIT*spaxel_area))**2, "Flux cal Variance", head)

    if input_parameters["debug"] == True:
        save_fits_cube(outFile_dark, noise_cube_dark, "dark noise variance", head)
        save_fits_cube(outFile_read_noise, noise_cube_read_noise**2, "read noise variance", head)
        save_fits_cube(outFile_ddetector_thermal, noise_cube_thermal, "detector thermal noise variance", head)

    # Calculate 5-sigma sensitivity
    sens_5sigma = 5.*np.median(noise_cube_total)/DIT*np.median(factor_calibration) # erg/s/cm2/um
    lcentral = np.median(output_lambs) # micron
    fnu = sens_5sigma*lcentral**2/(const.c.value*1e6) # erg/s/cm2/Hz
    sens_ABmag = -2.5*np.log10(fnu/3631./1e-23) # AB mag

    logging.info("Sensitivity point-source 5sigma = {:.2f} mag = {:.2e} erg/s/cm2/um at {:.3f} um".format(sens_ABmag, sens_5sigma, lcentral))

    noise_total = np.max(noise_cube_object) + np.median(n_observations*noise_cube_back + n_observations*noise_cube_dark + n_observations*noise_cube_thermal + n_observations*noise_cube_read_noise**2)

    fraction_noise_back = n_observations*np.median(noise_cube_back)/noise_total*100.
    fraction_noise_dark = n_observations*np.median(noise_cube_dark)/noise_total*100.
    fraction_noise_read = n_observations*np.median(noise_cube_read_noise)**2/noise_total*100.

    logging.info("Noise contributions: Background (sky+tel+instrument) = {:.2f} %. Dark current = {:.2f} %. Read noise = {:.2f} %".format(fraction_noise_back, fraction_noise_dark, fraction_noise_read))

    # CHANGELOG 11-01-2024: Removed the crosstalk application as we use visible detectors, not IR
    # # Save transmission with crosstalk
    # if grating != "V+R":
    #     if not input_parameters["mci"]:
    #         output_transmission = apply_crosstalk_1d(output_transmission, config_data["crosstalk"])

    np.savetxt(base_filename + "_total_tr.txt", np.c_[output_lambs, output_transmission], comments="#", header="\n".join([
        'TYPE: Total transmission.',
        'Wavelength [um], Transmission']))


    # Save PSF
    head_PSF = head.copy()
    head_PSF['CDELT1'] = spax_scale.psfscale
    head_PSF['CDELT2'] = spax_scale.psfscale
    head_PSF['LAMBDA'] = psf_lambda
    try:
        del head_PSF['CDELT3']
        del head_PSF['CRVAL3']
        del head_PSF['CRPIX3']
        del head_PSF['CTYPE3']
        del head_PSF['CUNIT3']
        del head_PSF['BUNIT']
        del head_PSF['SPECRES']
    except:
        pass

    save_fits_cube(base_filename + "_PSF_internal.fits", psf_internal, "PSF", head_PSF)

    # - Rebin and Reshape PSF to match the output cube size and spaxel size
    def rebin_psf(arr, new_shape):
        shape = (new_shape[0], arr.shape[0] // new_shape[0],
            new_shape[1], arr.shape[1] // new_shape[1])
        return arr.reshape(shape).sum(-1).sum(1)
    psf_info = config_data["spaxel_scale"][input_parameters["spaxel_scale"]]
    
    logging.info("psf_info: %s", psf_info) # CHANGELOG 09-01-2024: added for my troubleshooting

    # Calculate the offset needed to keep the PSF center
    # at the output image center after rebining
    psf_oversampling = int(round(min([psf_info.xscale, psf_info.yscale])/psf_info.psfscale))
    logging.info('psf_oversampling: %s', psf_oversampling)
    psf_spaxel_shape = psf_info.psfsize//psf_oversampling + 1

    tmp = np.zeros((psf_spaxel_shape*psf_oversampling, psf_spaxel_shape*psf_oversampling))
    psfcenter = psf_info.psfsize//2 + 1
    psfcenter_offset = psfcenter % psf_oversampling

    x0 = psf_oversampling//2 + 1 - psfcenter_offset
    tmp[x0:x0 + psf_info.psfsize, x0:x0 + psf_info.psfsize] = psf_internal[:, :]

    psf_spaxel_shape_x, psf_spaxel_shape_y = psf_spaxel_shape, psf_spaxel_shape

    psf_spaxel = rebin_psf(tmp, (psf_spaxel_shape, psf_spaxel_shape))

    if input_parameters["spaxel_scale"] == "30x60":
        # an extra rebin is needed for the y axis
        psf_spaxel_shape_y = psf_spaxel_shape_y//2
    elif input_parameters["spaxel_scale"] == "120x60":
        # an extra rebin is needed for the x axis
        psf_spaxel_shape_x = psf_spaxel_shape_x//2

    psf_spaxel = frebin2d(tmp, (psf_spaxel_shape_x, psf_spaxel_shape_y))
    psf_spaxel = psf_spaxel/np.sum(psf_spaxel)*np.sum(tmp) # Normalize PSF

    # center PSF on the output array
    center_x_output = (output_cube_spec.shape[2] - 1) // 2 - (output_cube_spec.shape[2] % 2 - 1)
    center_y_output = (output_cube_spec.shape[1] - 1) // 2 - (output_cube_spec.shape[1] % 2 - 1)
    center_x_psf = (psf_spaxel.shape[1] - 1) // 2 - (psf_spaxel.shape[1] % 2 - 1)
    center_y_psf = (psf_spaxel.shape[0] - 1) // 2 - (psf_spaxel.shape[0] % 2 - 1)

    # adjust y axis
    tmp = np.zeros((output_cube_spec.shape[1], psf_spaxel.shape[1]))
    if tmp.shape[0] > psf_spaxel.shape[0]:
        yi = center_y_output - center_y_psf
        tmp[yi:yi+psf_spaxel.shape[0], :] = psf_spaxel
    else:
        yi = center_y_psf - center_y_output
        tmp[:, :] = psf_spaxel[yi:yi+tmp.shape[0], :]
    psf_spaxel = tmp

    # adjust x axis
    tmp = np.zeros((output_cube_spec.shape[1], output_cube_spec.shape[2]))
    if tmp.shape[1] > psf_spaxel.shape[1]:
        xi = center_x_output - center_x_psf
        tmp[:, xi:xi+psf_spaxel.shape[1]] = psf_spaxel
    else:
        xi = center_x_psf - center_x_output
        tmp[:, :] = psf_spaxel[:, xi:xi+tmp.shape[1]]
    psf_spaxel = tmp

    head_PSF['CDELT1'] = spax_scale.xscale
    head_PSF['CDELT2'] = spax_scale.yscale
    save_fits_cube(base_filename + "_PSF.fits", psf_spaxel, "PSF", head_PSF)


    if msimlog.count_error == 0 and msimlog.count_warning == 0:
        logging.info('Simulation OK - ' + str(msimlog.count_error) + " errors and " + str(msimlog.count_warning) + " warnings")
    else:
        logging.warning('Simulation with problems - ' + str(msimlog.count_error) + " errors and " + str(msimlog.count_warning) + " warnings")


    logger.removeHandler(std)

    return 'Simulation complete'


def save_fits_cube(filename, data, typ, header):
    header['MSM_TYPE'] = typ
    fits.writeto(filename, data.astype(np.float32), header=header, overwrite=True, output_verify="silentfix")



class MSIMLoggingHandler(logging.Handler):
    def __init__(self):
        self.count_warning = 0
        self.count_error = 0
        logging.Handler.__init__(self)

    def emit(self, record):
        if record.levelno == logging.WARNING:
            self.count_warning += 1
        elif record.levelno >= logging.ERROR:
            self.count_error += 1


class MSIMFormatter(logging.Formatter):

    def __init__(self):
        logging.Formatter.__init__(self, "%(message)s")

    def format(self, record):
        if record.levelno > logging.INFO:
            self._style = logging._STYLES["%"][0]("** %(levelname)s ** %(message)s")
        else:
            self._style = logging._STYLES["%"][0]("%(message)s")

        return logging.Formatter.format(self, record)
