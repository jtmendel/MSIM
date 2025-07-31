'''
Calculates LSF, instrument background and transmission

CHANGELOG:

Version 0.0.0 (2023-10-27)
--------------------------
- Original HARMONI simulator code
Developers: Miguel Pereira Santaella, Laurence Routledge, Simon Zieleniewsk, Sarah Kendrew

Version 1.0.0 (2024-01-10)
--------------------------
- Changed the telescope area to be that of the telescope, not zero. (This was an issue from the HARMONI git!)
- Added logging info for troubleshooting.
- Changed logging info to use MAVIS rather than HARMONI.
- Updated FPRS and cryo temperature parameters to use MAVIS prefix.
- Removed MCI section of code.
Author: Eric Muller (eric.muller@anu.edu.au)

'''
import logging

import numpy as np
from scipy.interpolate import PchipInterpolator
import scipy.constants as sp
from astropy.convolution import Gaussian1DKernel
from astropy.io import fits

from src.config import *
from src.modules.misc_utils import path_setup
from src.modules.em_model import *


tppath = path_setup('../../' + config_data["data_dir"] + 'throughput/')
AreaTel = 0.


class InstrumentPart:
    """
    This is an updated InstrumentPart class specific to the MAVIS tabulated throughput format. 
    It expects that you provide a file describing the throughput of the part that you want to simulate
    and makes kludgy assumptions for emissivity (which shouldn't be a huge issue at MAVIS wavelengths).
    """
    edust = 0.5 # Grey Dust covering on some optics
    mindustfrac = 0.005 # 0.5% dust on optical surfaces - won't be perfectly clean

    def __init__(self, name, filename, temp, area, global_scaling=1.):
        
        self.name = name
        self.filename = filename
        self.temp = temp
        self.area = area
        self.global_scaling = global_scaling
        self.number = 0

    def set_number(self, number):
        self.number = number

    def calcThroughput(self, lamb, filename):
        # Read emissivity from filename
        try:
            l, tpt = np.loadtxt(os.path.join(tppath, filename), unpack=True, comments="#", delimiter=",")
        except:
            raise("Cannot load file: ", filename)

        # Interpolate throughput to output lambda grid
        tpt_interp = PchipInterpolator(l, tpt, extrapolate=True)
        
        return tpt_interp(lamb).clip(0,1)


    def calcThroughputAndEmission(self, lamb, DIT, output_file=""):
        throughput = self.calcThroughput(lamb, self.filename) * self.global_scaling
        
        #also get emissivity (kinda?)
        emissivity = (1. - throughput/self.global_scaling)
                
        #thermal nonsense
        emission = emissivity*blackbody(lamb, self.temp) #J/s/m2/lambda(um)/arcsec2
        
        #now as photons!
        emission_ph = emission/(sp.h*sp.c/(lamb*1.E-6))*DIT # photons/um/m2/arcsec2

        logging.debug("Instrument Part Model - {0}".format(self.name))
        logging.debug("global_scaling = {:5.3f}".format(self.global_scaling))
            
        logging.debug("lambda = {:7.4f} emissivity = {:6.3f} throughput = {:6.3f} emission_ph = {:.2e}".format(np.median(lamb), np.median(emissivity), np.median(throughput), np.median(emission_ph)))
        logging.debug("-------")
        
        if output_file is not None:
            plot_file = output_file + "_MAVIS_" + "{:02d}".format(self.number) + "_" + self.name.replace(" ", "_").lower()
                
            plt.clf()
            plt.plot(lamb, throughput)
            plt.xlabel(r"wavelength [$\mu$m]")
            plt.ylabel("Throughput " + self.name)
            plt.savefig(plot_file + "_tr.pdf")
            np.savetxt(plot_file + "_tr.txt", np.c_[lamb, throughput])

            plt.clf()
            plt.plot(lamb, emission_ph, label="Blackbody T = {:.1f} K".format(self.temp))
            plt.legend()
            plt.xlabel(r"wavelength [$\mu$m]")
            plt.ylabel("Emissivity " + self.name)
            plt.savefig(plot_file + "_em.pdf")
            np.savetxt(plot_file + "_em.txt", np.c_[lamb, emission_ph])
            
        
        return throughput, emission_ph


class Instrument:
    def __init__(self, name):
        self.name = name
        self.parts = []

    def addPart(self, part):
        self.parts.append(part)
        part.set_number(len(self.parts))

    def calcThroughputAndEmission(self, lamb, DIT, output_file=""):
        throughput = np.ones_like(lamb)
        emission = np.zeros_like(lamb)
    
        for part in self.parts:
            logging.info("- getting throughput for {0}".format(part.name))
            part_t, part_emi = part.calcThroughputAndEmission(lamb, DIT, output_file=output_file)
            
            # logging.info("part_t: %s, part_emi: %s", part_t, part_emi) # CHANGELOG 09-01-2024: added for my troubleshooting
  
            throughput *= part_t
            emission *= part_t
            emission = emission + part_emi

        # CHANGELOG 09-01-2024: Commented out the above, just print it regardless
        logging.info("Total MAVIS backround: lambda = {:7.4f}, throughput = {:6.3f}, emission = {:.2e} ph/um/m2/arcsec2/s"
                     .format(np.median(lamb), np.median(throughput), np.median(emission)/DIT)) # CHANGELOG 10-01-2024: Changed text from HARMONI to MAVIS


        return throughput, emission



def sim_instrument(input_parameters, cube, back_emission, transmission, ext_lambs, cube_lamb_mask, input_spec_res, debug_plots=False, output_file=""):
    ''' Simulates instrument effects
    Inputs:
        input_parameters: input dictionary
            exposure_time: Exposure time [s]
            grating: Spectral grating
            ao_mode: MCAO/NOAO/User defined PSF fits file
            telescope_temp: Telescope temperature [K]

        cube: Input datacube (RA, DEC, lambda)
        back_emission: Input background emission
        transmission: Input transmission
        ext_lambs: extended lambda array [um]
        cube_lamb_mask: mask array to get the lambs of the cube
        input_spec_res: Spectral resolution of the input cube [micron]
        debug_plots: Produce debug plots
        output_file: File name for debug plots

    Outputs:
        cube: cube including instrument effects
        back_emission: back_emission including telescope
        LSF_size: width of the LSF [A]
    '''
    
    # Get instrument transmission
    logging.info("Calculating MAVIS transmission and background") # CHANGELOG 10-01-2024: Changed text from HARMONI to MAVIS
    
    # CHANGELOG 11-01-2024: Changed to the MAVIS instrument
    mavis = Instrument("MAVIS")
    
    # Instrument model variables
    # -------------------------
    # Instrument temperatures (probably not super important given wavelength)
    TTel = input_parameters["telescope_temp"]
    TCool = 273.15 + config_data['MAVIS_FPRS_temp'] # fixed FPRS temp # CHANGELOG 29-12-2023: Changed to MAVIS prefix
    TCryo = config_data['MAVIS_cryo_temp'] # CHANGELOG 29-12-2023: Changed to MAVIS prefix
    TCryoMech = TCryo
    TTrap = TCool
    Touter_window = TTel - 0.2*(TTel - TCool)
    Tinner_window = TCool + 0.2*(TTel - TCool)
    AreaIns = (config_data['telescope']['diameter']*0.5)**2*np.pi    # Full aperture, including central obstruction
    global AreaTel
    AreaTel = config_data['telescope']['area']            # 8.2m with 1m central obscuration

    # Dust properties
    dustfrac = 0.01
    dustfrac = max(InstrumentPart.mindustfrac, dustfrac)    # Can make outer surfaces more dusty to represent aging

    # Cold trap properties
    ecoldtrap = 1.
    rwindow = 0.01    # 1% AR coating on each surface
    # -------------------------
    
    logging.debug("MAVIS model. TCool = {:d} K TCryo = {:d} K TCryoMech = {:d} K TTrap  = {:d} K Touter_window = {:d} K Tinner_window = {:d} K"
                    .format(*map(int, [TCool, TCryo, TCryoMech, TTrap, Touter_window, Tinner_window])))
    logging.debug("AreaIns = {:6.1f} m2 AreaTel = {:6.1f} m2".format(AreaIns, AreaTel))
    logging.debug("ecoldtrap = {:6.3f} rwindow = {:6.3f}".format(ecoldtrap, rwindow))
    logging.debug("-------")
    
    #included some kind of global scaling parameter to deal with budget vs. requirement sensitivity
    tpt_model = input_parameters['throughput_model']
    if tpt_model == 'budget':
        scale_factor = 1.
    elif tpt_model == 'requirement':
        scale_factor = 0.55

    #add some MAVIS parts!
    mavis.addPart(InstrumentPart("AO bench", "mavis_AOM_throughput_2025-03-14_spec.csv", TTel, AreaIns, global_scaling=scale_factor))

    #grating selection
    grating = input_parameters['grating']
    spaxel = input_parameters['spaxel_scale']
    if grating == 'HR-Blue':
        if spaxel == '25x25':
            mavis.addPart(InstrumentPart("Spectrograph + " + grating, "hrblue_grating_25mas_2025-03-06.csv", TTel, AreaIns, global_scaling=scale_factor))
        elif spaxel == '50x50':
            mavis.addPart(InstrumentPart("Spectrograph + " + grating, "hrblue_grating_50mas_2025-03-06.csv", TTel, AreaIns, global_scaling=scale_factor))
    elif grating == 'LR-Blue':
        if spaxel == '25x25':
            mavis.addPart(InstrumentPart("Spectrograph + " + grating, "lrblue_grating_25mas_2025-03-06.csv", TTel, AreaIns, global_scaling=scale_factor))
        elif spaxel == '50x50':
            mavis.addPart(InstrumentPart("Spectrograph + " + grating, "lrblue_grating_50mas_2025-03-06.csv", TTel, AreaIns, global_scaling=scale_factor))
    elif grating == 'HR-Red':
        if spaxel == '25x25':
            mavis.addPart(InstrumentPart("Spectrograph + " + grating, "hrred_grating_25mas_2025-03-06.csv", TTel, AreaIns, global_scaling=scale_factor))
        elif spaxel == '50x50':
            mavis.addPart(InstrumentPart("Spectrograph + " + grating, "hrred_grating_50mas_2025-03-06.csv", TTel, AreaIns, global_scaling=scale_factor))
    elif grating == 'LR-Red':
        if spaxel == '25x25':
            mavis.addPart(InstrumentPart("Spectrograph + " + grating, "lrred_grating_25mas_2025-03-06.csv", TTel, AreaIns, global_scaling=scale_factor))
        elif spaxel == '50x50':
            mavis.addPart(InstrumentPart("Spectrograph + " + grating, "lrred_grating_50mas_2025-03-06.csv", TTel, AreaIns, global_scaling=scale_factor))

    ################################ MAVIS frameworking ################################
    
    # Nasmyth flange -- spectrograph arms
    MAVIS_transmission, MAVIS_background = mavis.calcThroughputAndEmission(ext_lambs, input_parameters["exposure_time"], output_file=output_file)
    
    back_emission = back_emission*MAVIS_transmission
    transmission = transmission*MAVIS_transmission
    back_emission = back_emission + MAVIS_background
    
    # Add instrument emission/transmission to the input cube
    instrument_tr_cube = MAVIS_transmission[cube_lamb_mask]
    instrument_tr_cube.shape = (np.sum(cube_lamb_mask), 1, 1)
    cube *= instrument_tr_cube

    instrument_background_cube = MAVIS_background[cube_lamb_mask]
    instrument_background_cube.shape = (np.sum(cube_lamb_mask), 1, 1)
    cube += instrument_background_cube


    # - LSF
    logging.info("Convolve with LSF")
    # Assume Gaussian LSF
    bandws = config_data['gratings'][grating] # CHANGELOG 11-01-2024: Removed mci option
        
    new_res = (bandws.lmin + bandws.lmax)/(2.*bandws.R) # micron
    pix_size = (ext_lambs[1] - ext_lambs[0])
    if new_res > input_spec_res:
        new_res_pix = (new_res**2 - input_spec_res**2)**0.5/pix_size
    else:
        logging.warning("The output spectral resolution is higher than the input cube resolution. Assuming input resolution = 0 AA")
        new_res_pix = new_res/pix_size
        
    logging.info("Output resolution: {:.3f} AA".format(new_res*10000.))
    logging.info("Input resolution: {:.3f} AA".format(input_spec_res*10000.))
    logging.info("Effective LSF FWHM = {:.3f} AA".format(new_res_pix*pix_size*10000.))
    
    LSF_size = 0
    if new_res_pix > 1.: # avoid convolution with a kernel narrower than 1 pixel
        sigma_LSF_pix = new_res_pix/2.35482
    
        npix_LSF = int(sigma_LSF_pix*config_data['LSF_kernel_size'])
        # Ensure that the kernel has an odd number of channels
        if npix_LSF % 2 == 0:
            npix_LSF = npix_LSF + 1
            
        kernel_LSF = Gaussian1DKernel(stddev=sigma_LSF_pix, x_size=npix_LSF)
        z, y, x = cube.shape
        
        for py in range(y):
            for px in range(x):
                spectrum = np.copy(back_emission)
                spectrum[cube_lamb_mask] = cube[:, py, px]
                
                cube[:, py, px] = np.convolve(spectrum, kernel_LSF, mode="same")[cube_lamb_mask]
        
        
        back_emission = np.convolve(back_emission, kernel_LSF, mode="same")
        transmission = np.convolve(transmission, kernel_LSF, mode="same")
        
        LSF_size = npix_LSF*(ext_lambs[1] - ext_lambs[0])*10000. # AA
        logging.info("Range for the LSF convolution: {:.3f} AA".format(LSF_size))
    else:
        logging.warning("LSF convolution not performed because the effective LSF FWHM is < 1 pix")

    #this lives here just to keep things from breaking
    fpm_final = None
    
    return (cube, back_emission, transmission, fpm_final), LSF_size
