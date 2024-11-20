'''
Calculates detector QE, dark current, read noise

CHANGELOG:

Version 0.0.0 (2023-10-27)
--------------------------
- Original HARMONI simulator code
Developers: Miguel Pereira Santaella, Laurence Routledge, Simon Zieleniewsk, Sarah Kendrew

Version 1.0.0 (2024-01-10)
--------------------------
- Removed NIR saturation masking options and NIR detector performance options
Author: Eric Muller (eric.muller@anu.edu.au)

Version 1.0.1 (2024-)
- Changed the detector QE curve to the E2V290 CCD QE curve
- Changed the pixel size to 10 microns squared for the E2V290 CCD
- 
#TODO:


'''
import os
import logging

import numpy as np
import scipy.constants as sp
from scipy.interpolate import UnivariateSpline
from scipy.interpolate import PchipInterpolator as interp1d
from scipy import integrate
from astropy.convolution import Gaussian1DKernel
from astropy.io import fits

from src.modules.misc_utils import path_setup
from src.config import *
from src.modules.em_model import *

detpath = path_setup('../../' + config_data["data_dir"] + 'detectors/')

#Detector throughput curve generated just using wavelength array.
def detector_QE_curve(wavels, grating, debug_plots, output_file):
	'''Function that generates a detector QE curve.
	Current detector data taken from Detector Reference
	table.

	Inputs:
		wavels: array of wavelengths for datacube
		grating: grating choice to set detector

	Outputs:
		cube_det_qe: array of detector QE values for
			each wavelength in array
	'''
	# CHANGELOG 11-01-2024: For now, manually set to the E2V290 CCD QE curve.
	detector_QE_file = "E2V290_QE.txt"
		
	cube_det_qe, orig_qe_lambda, orig_qe = load_transmission_curve(wavels, detector_QE_file, debug_plots, [output_file, "det_qe"], "detector QE", scaling=1/100., full_curve=True)
	return cube_det_qe, orig_qe_lambda, orig_qe

def mask_saturated_pixels(cube, grating):
	''' Mask saturated pixels
	Inputs:
		cube: Input datacube (RA, DEC, lambda)
		grating: Spectral grating
	Outputs:
		cube: masked cube
		mask_pixels: saturated pixels mask
	'''
	logging.info("Masking saturated pixels")
	
	# CHANGELOG 03-01-2024: Just set it to the 'vis' values as we have only included these.
	limit = config_data["saturation"]["vis"]
	
	# CHANGELOG 03-01-2024: Commented this out, see above.
	# if grating == "V+R":
	# 	limit = config_data["saturation"]["vis"]
	# else:
	# 	limit = config_data["saturation"]["nir"]
		
	
	mask_pixels = cube > limit
	if np.sum(mask_pixels) > 0:
		logging.warning(str(np.sum(mask_pixels)) + " pixels are saturated in the observed cube")
		
	cube[mask_pixels] = 0.
	
	return cube, mask_pixels

def sim_detector(input_parameters, cube, back_emission, transmission, lambs, debug_plots=False, output_file=""):
	''' Simulates detector effects
	Inputs:
		input_parameters: input dictionary
			exposure_time: Exposure time [s]
			grating: Spectral grating

	
		cube: Input datacube (RA, DEC, lambda)
		back_emission: Input background emission
		transmission: Input transmission
		lambs: lambda array [um]
		grating: Spectral grating
		debug_plots: Produce debug plots
		output_file: File name for debug plots

	Outputs:
		cube: Cube including QE
		read_noise: read noise for the grating and DIT [e/pix]
		dark*DIT: dark current  [e/pix]
		ftot_electron*DIT: thermal background seen by the detector  [e/pix]
	'''
	
	DIT = input_parameters["exposure_time"]
	grating = input_parameters["grating"]
	
	detector_performance = config_data['detector']['avg']
	
	# Get QE curve
	logging.info("Calculating detector QE")
	if type(detector_performance['qe']) == str and detector_performance['qe'] == "file":
		qe_curve, orig_qe_lambda, orig_qe = detector_QE_curve(lambs, grating, debug_plots, output_file)
	else:
		f = interp1d(detector_performance['qe']["w"], detector_performance['qe']["qe"], extrapolate=True)
		qe_curve = f(lambs)
		
		if debug_plots:
			plot_file = [output_file, "det_qe"]
			plot_label = "detector QE"
			plt.clf()
			plt.plot(lambs, qe_curve)
			plt.xlabel(r"wavelength [$\mu$m]")
			plt.ylabel(plot_label)
			plt.savefig("_".join(plot_file) + "_tr.pdf")
			np.savetxt("_".join(plot_file) + "_tr.txt", np.c_[lambs, qe_curve])
			
		
		orig_qe_lambda = np.linspace(0.4, 2.6, 200)
		orig_qe = f(orig_qe_lambda)
		
		
	back_emission *= qe_curve
	transmission *= qe_curve
	
	qe_curve.shape = (len(lambs), 1, 1)
	cube *= qe_curve
 
	# CHANGELOG 03-01-2024: Just set it to the 'vis' values as we have only included these.
	read_noise = detector_performance["read_noise"]["vis"]
	dark = detector_performance["dark_current"]["vis"]

    # keeping this so that things don't break, but probably not correct. Thermal emission not really an issue anyways. 
	# Thermal emission seen by the detector
	# 10-micron pixels and F/2 camera... (independent of spaxel scale)
	pixel_area = 10E-6**2 # m2 # CHANGELOG 11-01-2024: Changed to 10um pixels for E2V290 CCD
	pixel_solid_cryo_mech = 0.2*(360/2./np.pi*3600.)**2 # arcsec**2
	pixel_solid_rad = 2.0*np.pi*(360/2./np.pi*3600.)**2 # arcsec**2
	
	TCryoMech = config_data["MAVIS_cryo_temp"] - 5.
	Trad = config_data["MAVIS_cryo_temp"]
	
	fcryomech = blackbody(orig_qe_lambda, TCryoMech)/(sp.h*sp.c/(orig_qe_lambda*1.E-6))*pixel_solid_cryo_mech # photons/s/um/m2
	frad = blackbody(orig_qe_lambda, Trad)/(sp.h*sp.c/(orig_qe_lambda*1.E-6))*pixel_solid_rad # photons/s/um/m2
	ftot = (fcryomech + frad)*pixel_area # photons/um/pixel
	ftot = ftot*orig_qe  # e/um/pixel
	ftot_electron = integrate.simpson(ftot, x=orig_qe_lambda) #e/pixel
	
	return cube, back_emission, transmission, read_noise, dark*DIT, ftot_electron*DIT
	

# Detector systematics code

def interp(data):
	"""
	Interpolation function to convert distribution into a callable
	function to extract read noise values from.

	Parameters:
			data - Data to draw from
	"""
	X = np.sort(data)
	x_range = np.linspace(0,1,len(X))
	func = UnivariateSpline(x_range, X, k=4, s=0)
	return func

def make_rn_dist(det_save_path):
	"""
	Draw random samples from distribution to create read
	noise values.

	Inputs:
		det_save_path: Directory to save rn distribution to
	"""
	if det_save_path == "None":
		det_save_path = detpath
	
	# First look to see if detectors already made
	if os.path.exists(det_save_path+'MAVIS_dets.fits'):
		logging.info('- found existing MAVIS detectors')
		if config_data['systematics']['force_new'] == False:
			logging.info('- using existing detectors')
			hdulist = fits.open(det_save_path+'MAVIS_dets.fits')
			rn_vals = hdulist[0].data
			return rn_vals
		else:
			logging.info('- overwriting exisiting detectors')

	logging.info('- no exisiting MAVIS detectors found')
	logging.info('- creating new detectors')
	
	rn_path = path_setup('../../' + config_data["data_dir"] + 'detectors/')
	rn_file = rn_path+config_data['systematics']['rn_file'] 
	if os.path.isfile(rn_file) is False:
		logging.error('There was an error finding rn_file. \
					  Check config.py to ensure filename is correct')
		return
	
    # TODO: change the readnoise file to that for the E2V290 CCD?
	hdu = fits.open(rn_file)
	rn_data = hdu[0].data

	side_len, N_det = config_data['side_length'], config_data['N_IR_det']
	N = (side_len**2) * N_det
	rands = np.random.random(N) 
	rands_sort = np.sort(rands)
	rn_vals = interp(rn_data.flatten())(rands_sort)
	
	np.random.shuffle(rn_vals)
	rn_vals = np.reshape(rn_vals,(N_det,side_len,side_len))

	#Saving detectors for future use
	fits.writeto(det_save_path+'HARMONI_dets.fits',rn_vals,overwrite=True)
	
	return rn_vals


def make_dets(rn_vals, DIT):
	"""
	Create simulated detector data for multiple read noise values.

	Inputs:
		rn_vals (list): List of read noise values for each detector.
		DIT (float): Detector integration time.

	Outputs:
		dets (ndarray): Array of simulated detector data for each detector.
		dets_head (Header): Header information of the last detector data.
	"""
	dets = []
	for i in range(8):
		ng_h4rg = ng.HXRGNoise()
		det_hdu = ng_h4rg.mknoise(o_file=None, rn_array=rn_vals[i], dit=DIT)
		det_data = det_hdu.data
		dets.append(det_data)
	dets = np.array(dets)
	dets_head = det_hdu.header
	return dets, dets_head


def add_detectors(cube, dets):
	''' Adds detectors to a datacube
	Inputs:
		cube: Input datacube (RA, DEC, lambda)
		dets: List of 2 detector arrays

	Outputs:
		new_datacube: Cube with detector systematics added
	'''
	sim_dets = np.copy(dets)
	cube_shape = cube.shape
	
	# allow for datacube to be smaller than the detectors
    # MAVIS 144x100, split 50+50
	#full_datacube = np.zeros((3700,152,204))
	full_datacube = np.zeros((8750,144,100))
	full_datacube[:cube_shape[0],:cube_shape[1],:cube_shape[2]] += cube

	# slice datacube up into 2 halves
	slice1 = full_datacube[:,:,:50]
	slice2 = full_datacube[:,:,50:]
	slices = [slice1, slice2]

	# add detector read noise onto datacube
	for i in range(len(sim_dets)):
		for j in range(slices[0].shape[1]):
			slitlet = slices[i][:,j,:]
			if i % 2 == 0:
				hoz_pos = np.int(np.round(32.5+(j*102)+(4.73*j)))
				if j % 2 == 0:
					vert_pos = 130 + (2 * j) - 67
				else:
					vert_pos = 130 + (2 * j) + 67
			else:
				hoz_pos = np.int(np.round(12.5+(j*102)+(4.73*j)))
				if j % 2 == 0:
					vert_pos = 204 - (2 * j) - 67
				else:
					vert_pos = 204 - (2 * j) + 67

			sim_dets[i][vert_pos:vert_pos+8750, hoz_pos:hoz_pos+50]+=slitlet

		# Convert detectors back into datacube		
	new_datacube_slices=[]
	for i in range(len(sim_dets)):
		new_datacube_slice = np.zeros(slices[i].shape)
		for j in range(new_datacube_slice.shape[1]):
			if i % 2 == 0:
				hoz_pos = np.int(np.round(32.5+(j*102)+(4.73*j)))
				if j % 2 == 0:
					vert_pos = 130 + (2 * j) - 67
				else:
					vert_pos = 130 + (2 * j) + 67
			else:
				hoz_pos = np.int(np.round(12.5+(j*102)+(4.73*j)))
				if j % 2 == 0:
					vert_pos = 100 - (2 * j) - 67 
				else:
					vert_pos = 100 - (2 * j) + 67
			new_datacube_slice[:,j,:] = sim_dets[i][vert_pos:vert_pos+8750, hoz_pos:hoz_pos+50]
		new_datacube_slices.append(new_datacube_slice)

	new_datacube = np.zeros((8750,144,100))

	new_datacube[:,:38,:102] = new_datacube_slices[0]
	new_datacube[:,:38,102:] = new_datacube_slices[1]
#	new_datacube[:,38:76,:102] = new_datacube_slices[2]
#	new_datacube[:,38:76,102:] = new_datacube_slices[3]
#	new_datacube[:,76:114,:102] = new_datacube_slices[4]
#	new_datacube[:,76:114,102:] = new_datacube_slices[5]
#	new_datacube[:,114:,:102] = new_datacube_slices[6]
#	new_datacube[:,114:,102:] = new_datacube_slices[7]

	# ensure new cube is same size as input cube
	new_datacube = new_datacube[:cube_shape[0],:cube_shape[1],:cube_shape[2]]
	return new_datacube
