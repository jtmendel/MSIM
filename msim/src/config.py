'''
File that stores hardwired data for use in MAVIS
simulation pipeline. Data stored in dictionary format
with keywords.

CHANGELOG:

Version 0.0.0 (2023-10-27)
--------------------------
- Original HARMONI simulator code
Developers: Miguel Pereira Santaella, Laurence Routledge, Simon Zieleniewsk, Sarah Kendrew

Version 1.0.0 (2024-01-10)
--------------------------
- Removed "best," "worst," and "contractual" detector config data of read noise, dark current, and QE
- Updated the read noise to 3 e/pix/hr as this is likely what MAVIS will operate at with its science frequency
- Updated the dark current to 3 e/pix/hr as this is likely what MAVIS will operate at 
- Updated the saturation values for the detectors. Removed the "nir" entry, updated the "vis" entry to 90,000e
- Updated the telescope mirror values to those for the VLT U4
- Renamed the FPRS_temp and cryo_temp variables to have the MAVIS_ prefix. Changed the cryo temp to the nominal detector operating temperature for MAVIS, 163K.
- Swapped over to the MAVIS grating values, from Table 1 of MAVIS-SPEC-REP-0002 V4. Removed nominal grating values.
- Changed the spaxel scales to those for MAVIS
- Removed MCI values from the code
Author: Eric Muller (eric.muller@anu.edu.au)

'''
import logging
import collections
import sys

GratingInfo = collections.namedtuple('GratingInfo', 'lmin, lmax, R')
SpaxelScaleInfo = collections.namedtuple('SpaxelScaleInfo', 'xscale, yscale, psfscale, psfsize')


config_data = {
    'detector':{
        "avg":{
            'read_noise': {"vis":2.5}, # e/pix 
            'dark_current': {"vis":4./3600}, 
            'qe':'file', 
        },
    },
    
    # CHANGELOG: 10-12-2023 - Removed nir, updated vis value. Using the average value from CCD spec sheet
    'saturation': {"vis":90000.}, # e
    
    'crosstalk': 0.02, # fraction of photons going to each of the 4 contiguous pixels
    #TODO: leave crosstalk as is for now, can be updated later.
    'side_length':4096, #TODO: figure out how to deal with this as we have different side lengths
    'N_IR_det':8, #TODO: remove, do we need a vis value?

    # JESSE, TREVOR: Any insight into these values and what they are?
    #Detector systematics parameters
#    'systematics': {"bias_offset":5000.,
#                    "bias_amp":500.,
#                    "force_new":False
#                   },
                        
    # JESSE, TREVOR: Any insight into these values?
    'spectral_sampling':{"output":2.3, "internal":4.}, # spectral sampling of the output cube and internal. Nyquist = 2
    'LSF_kernel_size':12., # LSF kernel size in sigma units
    
    # CHANGELOG: 19-12-2023 - Updated to the VLT values, from Jesse
    'telescope': {'diameter':7.83, 'area':48.1}, #diam [m], area [m^2], M2 acts as pupil stop so only ~8m effective diam., 
    
    # CHANGELOG: 10-12-2023 - Renamed from HARMONI prefix
    'MAVIS_FPRS_temp': +2., # C #TODO: find where this is referenced, as we have changed from HARMONI prefix to MAVIS prefix for variable
    'MAVIS_cryo_temp': 163., # K # CHANGELOG 10-01-2024: changed to the nominal operating termperature of MAVIS
    
    'data_dir':"sim_data/",
    
    # CHANGELOG: 19-12-2023 - Changed over to MAVIS grating values, from MAVIS-SPEC-REP-0002 V4, Table 1. Note we leave them in um as the HARMONI values
    'gratings': {
            'LR-Blue':GratingInfo(0.370, 0.720, 4000.),
            'LR-Red':GratingInfo(0.510, 0.950, 4000.),
            'HR-Blue':GratingInfo(0.425, 0.550, 12000.),
            'HR-Red':GratingInfo(0.630, 0.880, 9000.)
            },

    #TODO: check what the values of psfscale and psfsize should be and what they do in the code
    # CHANGELOG: 19-12-2023 - Switched to the 25 and 50 mas spaxel sizes. Need to figure out what the values for the psfscale and psfsize. Arbitrary values of 1.0 and 1000.0 have been used for now.
    'spaxel_scale': { '25x25':SpaxelScaleInfo(25., 25., 5., 1024),
                  '50x50':SpaxelScaleInfo(50., 50., 5., 1024)
                  },
    
    
    #FWHM of Instrument PSF depending on output spaxel scale in mas
    #Factors taken into account:
    #design image quality, manufacturing and assembly tolerances, vibration, flexure, diff refraction,
    'dynamic_instrument_psf': 5.5,
    # CHANGELOG 29-12-2023: Changed to the MAVIS pixel scales, used arbitrary values for now.
    #TODO: check what values we should be using
    'static_instrument_psf': {
            '25x25': 5.,
            '50x50':5.,
        },
    #TODO: above will come from optical model
    # CHANGELOG: 10-12-2023 - Removed minimum compliant instrument values ('mci_dynamic_instrument_psf' and 'mci_static_instrument_psf')

    # Each PSD file containts 1 seeing  [0.43] and 1 zenith angle [25] #TODO: ask jesse if we need these?
    'PSD_file':{"MCAO":"psf_mcao_msim_8LGS.fits.fz"},
    'PSD_params':{'air_masses':[1.1],
        'seeings':[0.80]} # arcseconds
    #'PSD_params':{'air_masses':[1.15, 1.3, 1.5, 2.0],
    #    'seeings':[0.80, 0.57, 0.64, 0.72, 1.04]} # arcseconds

}

class MSIMError(Exception):
    pass
    def __init__(self, message):
        logging.error(message)
        sys.exit()

