'''
Calculates LSF, instrument background and transmission
'''
import logging

import numpy as np
from scipy.interpolate import interp1d, interp2d
import scipy.constants as sp
from astropy.convolution import Gaussian1DKernel
from astropy.io import fits

from src.config import *
from src.modules.misc_utils import path_setup
from src.modules.em_model import *


tppath = path_setup('../../' + config_data["data_dir"] + 'throughput/')
hc_path = path_setup('../../' + config_data["data_dir"] + 'HC/')

class InstrumentPart:
	"""
	Represents a part of an instrument.

	Attributes:
		substrate (str): The substrate material.
		mirror (str): The mirror material.
		edust (float): The amount of grey dust covering on some optics.
		mindustfrac (float): The percentage of dust on optical surfaces.

	Methods:
		__init__(self, name, temp, area, n_mirrors=0, n_lenses=0, dust_lens=0., dust_mirror=mindustfrac, global_scaling=1., emis_scaling=1., emis_mirror=mirror, emis_lens=substrate, emis_dust=edust): Initializes the InstrumentPart object.
		set_number(self, number): Sets the number of the instrument part.
		calcEmissivity(self, lamb, filename, scaling, dust, n): Calculates the emissivity of the instrument part.
		calcThroughputAndEmission(self, lamb, DIT, output_file=""): Calculates the throughput and emission of the instrument part.
	"""
	substrate = "Suprasil3001_50mm_Emissivity.txt"
	mirror = "QuantumFS500_Emissivity.txt"
	edust = 0.5
	mindustfrac = 0.005

	def __init__(self, name, temp, area, n_mirrors=0, n_lenses=0, dust_lens=0., dust_mirror=mindustfrac, global_scaling=1., emis_scaling=1., emis_mirror=mirror, emis_lens=substrate, emis_dust=edust):
		"""
		Initialize the InstrumentPart object.

		Args:
			name (str): The name of the instrument.
			temp (float): The temperature of the instrument.
			area (float): The area of the instrument.
			n_mirrors (int, optional): The number of mirrors. Defaults to 0.
			n_lenses (int, optional): The number of lenses. Defaults to 0.
			dust_lens (float, optional): The amount of dust on lenses. Defaults to 0.
			dust_mirror (float, optional): The amount of dust on mirrors. Defaults to mindustfrac.
			global_scaling (float, optional): The global scaling factor. Defaults to 1.
			emis_scaling (float, optional): The emissivity scaling factor. Defaults to 1.
			emis_mirror (float, optional): The emissivity of mirrors. Defaults to mirror.
			emis_lens (float, optional): The emissivity of lenses. Defaults to substrate.
			emis_dust (float, optional): The emissivity of dust. Defaults to edust.
		"""
		if dust_lens != 0:
			assert(n_lenses != 0)
		
		self.name = name
		self.temp = temp
		self.area = area
		self.n_mirrors = n_mirrors
		self.n_lenses = n_lenses
		self.dust_lens = dust_lens
		self.dust_mirror = dust_mirror
		self.global_scaling = global_scaling
		self.emis_scaling = emis_scaling
		self.emis_mirror = emis_mirror
		self.emis_lens = emis_lens
		self.emis_dust = emis_dust
		self.number = 0

	def set_number(self, number):
		"""
		Sets the number of the instrument part.

		Args:
			number (int): The number of the instrument part.
		"""
		self.number = number

	def calcEmissivity(self, lamb, filename, scaling, dust, n):
		"""
		Calculates the emissivity of the instrument part.

		Args:
			lamb (float): The wavelength.
			filename (str or float): The filename or number representing the emissivity.
			scaling (float): The scaling factor.
			dust (float): The amount of dust.
			n (int): The number of elements.

		Returns:
			float: The calculated emissivity.
		"""
		# Code implementation
		#TODO: implement this function? Wasn't here as of version 311

	def calcThroughputAndEmission(self, lamb, DIT, output_file=""):
		"""
		Calculates the throughput and emission of the instrument part.

		Args:
			lamb (float): The wavelength.
			DIT (float): The DIT value.
			output_file (str, optional): The output file name. Defaults to "".

		Returns:
			tuple: A tuple containing the throughput and emission values.
		"""
		# Code implementation
		#TODO: implement this function? Wasn't here as of version 311		
	
class Instrument:
	"""Represents an instrument.

	Args:
		name (str): The name of the instrument.

	Attributes:
		name (str): The name of the instrument.
		parts (list): A list of parts in the instrument.

	"""

	def __init__(self, name):
		self.name = name
		self.parts = []

	def addPart(self, part):
		"""Adds a part to the instrument.

		Args:
			part (object): The part to be added.

		"""
		self.parts.append(part)
		part.set_number(len(self.parts))

	def calcThroughputAndEmission(self, lamb, DIT, output_file=""):
		"""Calculates the throughput and emission of the instrument.

		Args:
			lamb (numpy.ndarray): The wavelength array.
			DIT (float): The detector integration time.
			output_file (str, optional): The output file path. Defaults to "".

		Returns:
			tuple: A tuple containing the calculated throughput and emission arrays.

		"""
		throughput = np.ones_like(lamb)
		emission = np.zeros_like(lamb)

		for part in self.parts:
			part_t, part_emi = part.calcThroughputAndEmission(lamb, DIT, output_file=output_file)

			throughput *= part_t
			emission *= part_t
			emission = emission + part_emi

		if output_file is not None:
			logging.info("Total HARMONI backround: lambda = {:7.4f} throughput = {:6.3f} emission = {:.2e} ph/um/m2/arcsec2/s".format(np.median(lamb), np.median(throughput), np.median(emission)/DIT))

		return throughput, emission



def sim_instrument(input_parameters, cube, back_emission, transmission, ext_lambs, cube_lamb_mask, input_spec_res, debug_plots=False, output_file=""):
	''' Simulates instrument effects
	Inputs:
		input_parameters: input dictionary
			exposure_time: Exposure time [s]
			grating: Spectral grating
			ao_mode: LTAO/SCAO/NOAO/AIRY/User defined PSF fits file
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
	logging.info("Calculating HARMONI transmission and background")
	
	
	harmoni = Instrument("HARMONI")
	
	# Instrument model variables
	# -------------------------
	# Instrument temperatures
	TTel = input_parameters["telescope_temp"]
	#TCool = TTel - config_data['HARMONI_FPRS_diff_temp']
	TCool = 273.15 + config_data['HARMONI_FPRS_temp'] # fixed FPRS temp
	TCryo = config_data['HARMONI_cryo_temp']
	TCryoMech = TCryo
	TTrap = TCool
	Touter_window = TTel - 0.2*(TTel - TCool)
	Tinner_window = TCool + 0.2*(TTel - TCool)
	AreaIns = (config_data['telescope']['diameter']*0.5)**2*np.pi	# Full 37m2 aperture, including central obstruction -- this what we see from a thermal point of view after cold stop
	AreaTel = config_data['telescope']['area']			# 37m with 11m central obscuration -- this is what we see before the cold stop

	# Dust properties
	dustfrac = 0.01
	dustfrac = max(InstrumentPart.mindustfrac, dustfrac)	# Can make outer surfaces more dusty to represent aging

	# Cold trap properties
	ecoldtrap = 1.
	rwindow = 0.01	# 1% AR coating on each surface
	# -------------------------
	
	logging.debug("HARMONI model. TCool = {:d} K TCryo = {:d} K TCryoMech = {:d} K TTrap  = {:d} K Touter_window = {:d} K Tinner_window = {:d} K".format(*map(int, [TCool, TCryo, TCryoMech, TTrap, Touter_window, Tinner_window])))
	logging.debug("AreaIns = {:6.1f} m2 AreaTel = {:6.1f} m2".format(AreaIns, AreaTel))
	logging.debug("edust = {:6.3f} dustfrac = {:6.3f} mindustfrac = {:6.3f}".format(InstrumentPart.edust, dustfrac, InstrumentPart.mindustfrac))
	logging.debug("ecoldtrap = {:6.3f} rwindow = {:6.3f}".format(ecoldtrap, rwindow))
	logging.debug("-------")
	
	
	# AO dichroic if present
	aoMode = input_parameters["ao_mode"]
	if aoMode == "LTAO":
		harmoni.addPart(InstrumentPart("LTAO dichroic", TTel, AreaIns, n_lenses=1, emis_lens="LTAO_0.6_dichroic.txt", dust_lens=2.*dustfrac))
		harmoni.addPart(InstrumentPart("AO cold trap", TTrap, AreaIns, n_mirrors=1, emis_mirror=0., dust_mirror=0.03, emis_dust=ecoldtrap))
	
	harmoni.addPart(InstrumentPart("Outer window", Touter_window, AreaIns, n_lenses=1, emis_scaling=0.5, dust_lens=dustfrac + InstrumentPart.mindustfrac))
	harmoni.addPart(InstrumentPart("Inner window", Tinner_window, AreaIns, n_lenses=1, emis_scaling=0.5, dust_lens=2.*InstrumentPart.mindustfrac))

	harmoni.addPart(InstrumentPart("Window reflected", TTrap, AreaIns, n_mirrors=1, emis_mirror=0., dust_mirror=2.*0.8*2.0*rwindow, emis_dust=ecoldtrap))

	low_dust_iso6 = 1
	# FPRS
	if low_dust_iso6:
		harmoni.addPart(InstrumentPart("FPRS", TCool, AreaTel, n_mirrors=4))
	else:
		harmoni.addPart(InstrumentPart("FPRS", TCool, AreaTel, n_mirrors=3, dust_mirror=0.032))
		harmoni.addPart(InstrumentPart("FPRS", TCool, AreaTel, n_mirrors=1, dust_mirror=0.008))
	
	if aoMode in ["SCAO", "HCAO"]:
		harmoni.addPart(InstrumentPart("SCAO dichroic", TCool, AreaIns, n_lenses=1, emis_lens="SCAO_0.8_dichroic.txt", dust_lens=2.*dustfrac))
	
	if low_dust_iso6:
		harmoni.addPart(InstrumentPart("Cryo window", TCool, AreaTel, n_lenses=1, emis_scaling=0.4, dust_lens=InstrumentPart.mindustfrac))
	else:
		harmoni.addPart(InstrumentPart("Cryo window", TCool, AreaTel, n_lenses=1, emis_scaling=0.4, dust_lens=0.125))
		
	harmoni.addPart(InstrumentPart("Cryo window inner dust", TCryo+50., AreaIns, n_mirrors=1, emis_mirror=0., dust_mirror=InstrumentPart.mindustfrac))
	harmoni.addPart(InstrumentPart("Cryo window cold trap", TCryo+50., AreaIns, n_mirrors=1, emis_mirror=0., dust_mirror=2.0*rwindow, emis_dust=ecoldtrap))

	# Cryostat
	harmoni.addPart(InstrumentPart("Pre-optics+IFU+Spectrograph", TCryoMech, AreaIns, n_lenses=11, n_mirrors=23))

	# Grating
	grating = input_parameters["grating"]
	harmoni.addPart(InstrumentPart("Grating " + grating, TCryoMech, AreaIns, n_mirrors=1, emis_mirror=grating + "_grating.txt", dust_mirror=0))
	
	HARMONI_transmission, HARMONI_background = harmoni.calcThroughputAndEmission(ext_lambs, input_parameters["exposure_time"], output_file=output_file)
	
	if input_parameters["mci"]:
		logging.info("Using minimum compliant instrument background and throughput")
		
		bandws = config_data['gratings'][grating]
		mci_lamb = np.arange(bandws.lmin, bandws.lmax, (bandws.lmax+bandws.lmin)*0.5/bandws.R)
		mci_HARMONI_transmission, _ = harmoni.calcThroughputAndEmission(mci_lamb, input_parameters["exposure_time"], output_file=None)
		
		scaling_transmission = 1./np.mean(mci_HARMONI_transmission)*0.26
		HARMONI_transmission = np.interp(ext_lambs, mci_lamb, mci_HARMONI_transmission*scaling_transmission)
		
		#scaling_background = 2.031e3/(np.median(HARMONI_background)/input_parameters["exposure_time"])*np.median(HARMONI_transmission)/0.289
		scaling_background = 1.0754 # scaled background to match specifictaion at 2.2um at 9C
		HARMONI_background = scaling_background*HARMONI_background
		

		logging.info("Total MCI HARMONI backround: lambda = {:7.4f} throughput = {:6.3f} emission = {:.3e} ph/um/m2/arcsec2/s".format(np.median(ext_lambs), np.median(HARMONI_transmission), np.median(HARMONI_background)/input_parameters["exposure_time"]))

		plot_file = output_file + "_HARMONI_mci"
		
		plt.clf()
		plt.plot(ext_lambs, HARMONI_transmission)
		plt.xlabel(r"wavelength [$\mu$m]")
		plt.ylabel("Throughput mci")
		plt.savefig(plot_file + "_tr.pdf")
		np.savetxt(plot_file + "_tr.txt", np.c_[ext_lambs, HARMONI_transmission])

		plt.clf()
		plt.plot(ext_lambs, HARMONI_background, label="HARMONI mci")
		plt.legend()
		plt.xlabel(r"wavelength [$\mu$m]")
		plt.ylabel("Emissivity mci")
		plt.savefig(plot_file + "_em.pdf")
		np.savetxt(plot_file + "_em.txt", np.c_[ext_lambs, HARMONI_background])
		
		
		
		
		

	back_emission = back_emission*HARMONI_transmission
	transmission = transmission*HARMONI_transmission
	back_emission = back_emission + HARMONI_background
	
	# Add instrument emission/transmission to the input cube
	
	instrument_tr_cube = HARMONI_transmission[cube_lamb_mask]
	instrument_tr_cube.shape = (np.sum(cube_lamb_mask), 1, 1)
	cube *= instrument_tr_cube

	instrument_background_cube = HARMONI_background[cube_lamb_mask]
	instrument_background_cube.shape = (np.sum(cube_lamb_mask), 1, 1)
	cube += instrument_background_cube


	# - LSF
	logging.info("Convolve with LSF")
	# Assume Gaussian LSF
	if input_parameters["mci"]:
		bandws = config_data['gratings_nominal'][grating]
	else:
		bandws = config_data['gratings'][grating]
		
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


	# Apply high-constrast focal plane mask
	if aoMode == "HCAO":
		logging.info("Apply HC focal plane mask " + input_parameters["hc_fp_mask"])
		fpm = fits.getdata(os.path.join(hc_path, input_parameters["hc_fp_mask"] + ".fits.gz"), 0, memmap=True) # 0.39 mas sampling
		fpm_sampling = 0.39 # mas
		y, x = fpm.shape
		
		mask_xsize = x*fpm_sampling
		mask_ysize = y*fpm_sampling

		spax = input_parameters["spaxel_scale"]
		pix_size = config_data["spaxel_scale"][spax].psfscale
		cube_xsize = cube.shape[2]*pix_size
		cube_ysize = cube.shape[1]*pix_size
		
		xgrid_in = np.linspace(-abs(mask_xsize)*0.5, abs(mask_xsize)*0.5, x)
		ygrid_in = np.linspace(-abs(mask_ysize)*0.5, abs(mask_ysize)*0.5, y)

		xgrid_out = np.arange(-abs(cube_xsize)*0.5, abs(cube_xsize)*0.5, abs(pix_size))
		ygrid_out = np.arange(-abs(cube_ysize)*0.5, abs(cube_ysize)*0.5, abs(pix_size))

		fpm_interp = interp2d(xgrid_in, ygrid_in, fpm, kind='linear')
		fpm_final = fpm_interp(xgrid_out, ygrid_out)
		
		for i in range(cube.shape[0]):
			cube[i,:,:] *= fpm_final
	
	else:
		fpm_final = None

	
	return (cube, back_emission, transmission, fpm_final), LSF_size
