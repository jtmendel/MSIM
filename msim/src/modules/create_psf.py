'''
Module to create PSF from PSD

CHANGELOG:

Version 0.0.0 (2023-10-27)
--------------------------
- Original HARMONI simulator code
Developers: Miguel Pereira Santaella, Laurence Routledge, Simon Zieleniewsk, Sarah Kendrew

'''
import os
import logging

import numpy as np
from scipy.interpolate import PchipInterpolator
import scipy.ndimage
from astropy.io import fits
from scipy.interpolate import interp2d
from scipy.signal import fftconvolve

try:
    from ..config import *
except:
    from config import *

try:
    from src.modules.misc_utils import path_setup
    from src.modules.rebin import *
except:
    from modules.misc_utils import path_setup
    from modules.rebin import *

psf_path = path_setup('../../' + config_data["data_dir"] + 'PSF/')


# Thierry FUSCO 18/12/17 17:42
# SIMUL_PSF_DataPackage.zip
#Programme pour prendre en compte la multi-analyse les geometries
#d'etoiles et la postion de la galaxie
def psd_to_psf(psd, pup, D, phase_static = None, samp = None, fov = None, lamb = 2.2*1.e-6, jitter=0.):
    """
    Computes the Point Spread Function (PSF) from a residual phase Power Spectral Density (PSD) and a pupil shape.

    Args:
        psd (ndarray): 2D array with PSD values (in nm^2 per frequency at the PSF wavelength).
        pup (ndarray): 2D array representing the pupil.
        D (float): Pupil diameter.
        phase_static (ndarray, optional): Static phase in nm. Defaults to None.
        samp (float, optional): Final PSF sampling (number of pixels in the diffraction). Min = 2! Defaults to None.
        fov (float, optional): PSF Field of View (in arcsec). Defaults to None.
        lamb (float, optional): PSF wavelength in meters. Defaults to 2.2*1.e-6.
        jitter (float, optional): Gaussian jitter in pixels. Defaults to 0.

    Returns:
        ndarray: The computed PSF.

    Raises:
        MSIMError: If the PSD horizon is not at least two times larger than the pupil diameter.
        MSIMError: If overFoV != 1 (not fully tested).

    """ 
 
    dim = psd.shape[0]
    npup = pup.shape[0]
    
    sampnum = float(dim)/npup # numerical sampling related to PSD vs pup dimension
    L = D*sampnum # Physical size of the PSD
    
    if dim < 2*npup:
        raise MSIMError("the PSD horizon must at least two time larger than the pupil diameter")
    
    convnm = (2*np.pi/(lamb*1e9)) # nm to rad

    # from PSD to structure function
    Bg = np.fft.fft2(np.fft.fftshift(psd)*convnm**2)/L**2
    
    # creation of the structure function
    Dphi = np.real(2*(Bg[0, 0]-Bg))
    Dphi = np.fft.fftshift(Dphi)

    if samp is not None:
        sampin = samp
    else:
        sampin = sampnum
    
    if float(sampin) <= sampnum:
        dimnum = float(int(dim*(sampin/sampnum)/2)*2) # even dimension of the num psd
        sampout = dimnum / npup # real sampling
        Dphi2 = Dphi[int(dim/2-sampout*npup/2):int(dim/2+sampout*npup/2),int(dim/2-sampout*npup/2):int(dim/2+sampout*npup/2)]
        #print 'input sampling = ', str(sampin) , ' ---  output sampling = ', str(sampout),' --- max num sampling = ', str(sampnum)
    else:
        dimnum = float(int(dim*(sampin/sampnum)/2)*2) # even dimension of the num psd
        sampout = dimnum/npup
        Dphi2 = np.zeros((int(dimnum), int(dimnum)))+(Dphi[0,0]+Dphi[int(dim)-1, int(dim)-1]+Dphi[0, int(dim)-1]+Dphi[int(dim-1), 0])/4.
        Dphi2[int(dimnum/2-dim/2):int(dimnum/2+dim/2),int(dimnum/2-dim/2):int(dimnum/2+dim/2)] = Dphi
        #print 'WARNING : Samplig > Dim DSP / Dim pup => extrapolation !!! We rceommmend to increase the PSD size'
        #print 'input sampling = ', str(sampin) , ' ---  output sampling = ', str(sampout),' --- max num sampling = ', str(sampnum)

    
    # increasing the FoV PSF means oversampling the pupil
    FoVnum = 1./2.*(lamb/(sampnum*D))*dim/(4.85*1.e-6)
    if fov is None:
        fov = FoVnum
    
    overFoV = fov/FoVnum

    if overFoV != 1:
        raise MSIMError("overFoV != 1 not fully tested")
        dimover = float(int(dimnum*overFoV/2)*2)
        npupover =  float(int(npup*overFoV/2)*2)
        xxover  = np.arange(dimover)/dimover*dimnum
        xxpupover  = np.arange(npupover)/npupover*npup
        
        fDphi2 = interp2d(np.arange(Dphi2.shape[0]), np.arange(Dphi2.shape[0]), Dphi2, kind='cubic')
        Dphi2 = fDphi2(xxover, xxover)
        Dphi2[Dphi2 < 0] = 0.
        
        fpup = interp2d(np.arange(pup.shape[0]), np.arange(pup.shape[0]), pup, kind='cubic')
        pupover = fpup(xxpupover, xxpupover)
        pupover[pupover < 0] = 0.
    else:
        dimover = dimnum
        npupover = npup
        pupover = pup
    
    #print "dimover = ", dimover
    
    if phase_static is not None:
        npups    = phase_static.shape[0]
        if npups != npup:
            raise MSIMError("pup and static phase must have the same number of pixels")

        if overFoV != 1:
            fphase_static = interp2d(np.arange(phase_static.shape[0]), np.arange(phase_static.shape[0]), phase_static, kind='cubic')
            phase_static_o = fphase_static(xxpupover, xxpupover)
            #phase_static_o[phase_static_o < 0] = 0.
            
        else:
            phase_static_o = phase_static

    #print 'input FoV = ', str(fov) , ' ---  output FoV = ', str(FoVnum*float(dimover)/dimnum), ' ---  Num FoV = ', str(FoVnum)

    #if fov > 2*FoVnum:
        #print 'Warning : Potential alisiang issue .. I recommend to create initial PSD and pupil with a larger numbert of pixel'


    # creation of a diff limited OTF (pupil autocorrelation)
    tab = np.zeros((int(dimover), int(dimover)), dtype=complex)
    if phase_static is None:
        tab[0:pupover.shape[0], 0:pupover.shape[1]] = pupover
    else:
        tab[0:pupover.shape[0], 0:pupover.shape[1]] = pupover*np.exp(1j*phase_static_o*2*np.pi/lamb)
    

    dlFTO = np.real(np.fft.ifft2(np.abs(np.fft.fft2(tab))**2))
    dlFTO = np.fft.fftshift(np.abs(dlFTO)/np.sum(pup))
    
    # creation of AO OTF
    aoFTO     = np.exp(-Dphi2/2.)
    
    # Computation of final OTF
    sysFTO = aoFTO*dlFTO
    sysFTO = np.fft.fftshift(sysFTO)

    # add Gaussian jitter
    if np.sum(jitter) > 0.:
        sigmax = 1./(2.*np.pi*jitter[0])*sysFTO.shape[0]
        sigmay = 1./(2.*np.pi*jitter[1])*sysFTO.shape[1]
        
        
        Gauss2D = lambda x, y: 1./(2.*np.pi*sigmax*sigmay)*np.exp(-0.5*((x/sigmax)**2 + (y/sigmay)**2))
        xgrid = np.linspace(1, sysFTO.shape[0], sysFTO.shape[0]) - sysFTO.shape[0]*0.5 - 0.5
        ygrid = np.linspace(1, sysFTO.shape[1], sysFTO.shape[1]) - sysFTO.shape[1]*0.5 - 0.5
        xx, yy = np.meshgrid(xgrid, ygrid)
        kernel = np.fft.fftshift(Gauss2D(xx, yy))
        sysFTO = sysFTO*kernel

    # simulate the rotation of the PSF
    #if rotation_angle is not None:
        #nsteps = 15
        #angles = np.linspace(0, rotation_angle, nsteps)
        
        #tmpFTO = np.fft.fftshift(sysFTO)
        
        #sysFTO_tmp = np.zeros_like(tmpFTO)
        #for a in angles:
            #sysFTO_tmp += scipy.ndimage.rotate(tmpFTO, a, reshape=False)

        #sysFTO_tmp /= nsteps
        #sysFTO = np.fft.fftshift(sysFTO_tmp)
        
    
    # Computation of final PSF
    sysPSF = np.real(np.fft.fftshift((np.fft.fft2(sysFTO))))
    sysPSF = sysPSF/np.sum(sysPSF) #normalisation to 1

    return sysPSF

psfscale = None
fov = None
AO_mode = None

rotation_angle = None

# AO variables
pup = None
stats = None
psd = None
xgrid_out = None
ygrid_out = None
jitter = None
diameter = None

# no AO variables
zenith_seeing = None
air_mass = None

# user-defined PSF
user_psf = None

def set_jitter(_jitter):
    global jitter
    jitter = _jitter

def define_psf(input_parameters, _jitter, _fov, _psfscale, rotation=None):
    '''
    Define parameters used for the PSF generation

    Args:
        input_parameters:
            ao_mode: AO mode LTAO, SCAO, Airy, user
            ao_star_hmag: H magnitude of the LTAO AO star
            ao_star_distance: Distance from HARMONI FoV to LTAO AO star
            air_mass: Air mass of the observation
            zenith_seeing: Atmospheric seeing FWHM [arcsec]
            user_defined_psf: FITS file with the user defined PSF
        
        _jitter: Residual usre defined jitter and instrument PSF effect [mas]
        fov: number of pixels of the PSF
        psfscale: pixel size for the PSF [mas]
        
    Returns:
        None
    '''
    global pup, stats, psd, xgrid_out, ygrid_out, jitter, psfscale, fov, diameter, AO_mode, rotation_angle
    global zenith_seeing, air_mass
    global user_psf
    
    AO_mode = input_parameters["ao_mode"].upper()
    rotation_angle = rotation
    
    fov = _fov
    psfscale = _psfscale
    xgrid_out = (np.linspace(0, fov-1, fov) - fov*0.5)*psfscale
    ygrid_out = (np.linspace(0, fov-1, fov) - fov*0.5)*psfscale

    zenith_seeing = input_parameters["zenith_seeing"]
    air_mass = input_parameters["air_mass"]
 
#    if AO_mode in ["AIRY", "MCAO"]:
#        jitter = _jitter
#        diameter = config_data["telescope"]["diameter"]
#        logging.info("define AO PSF - " + AO_mode)
#        
#        # not using any static phase map
#        stats = None
#       
#        #VLT pupil image
#        pup = fits.getdata(os.path.join(psf_path, "ut4pupil640p.fits"))
#       
#        if AO_mode == "AIRY":
#            psd = np.zeros((2048, 2048))
#        else:
#            logging.info("Using PSD file: " + config_data["PSD_file"][AO_mode])
#            psd = fits.getdata(os.path.join(psf_path, config_data["PSD_file"][AO_mode]))
#
#            jitter_PSD = 2.
#
#            logging.info("Total AO jitter = {:.2f} mas".format(jitter_PSD))
#        
#            # combine PSD jitter and instrument and extra user defined jitter
#            jitter = (jitter**2 + jitter_PSD**2)**0.5
#            logging.info("Total PSF jitter = {0:.2f}x{1:.2f} mas".format(*jitter))
#
   
    if AO_mode in ["AIRY"]:
        jitter = _jitter
        diameter = config_data["telescope"]["diameter"]
        logging.info("define AO PSF - " + AO_mode)
        
        # not using any static phase map
        stats = None
       
        #VLT pupil image
        pup = fits.getdata(os.path.join(psf_path, "ut4pupil640p.fits"))
        psd = np.zeros((2048, 2048))

    elif AO_mode == "NOAO":
        logging.info("define noAO Gaussian PSF")
        
    elif AO_mode == "MCAO":
        #JTM modified to work with PSF cube (wave, x, y) to deal with wavelength
        #dependence. This means that user_psf is now an interpolator object, not an image.

        #need to parse this for dependence on the scale changer
        spaxel = input_parameters['spaxel_scale']
        if spaxel == '25x25':
            psffile = 'psf_mcao_msim_25mas.fits.fz'
        elif spaxel == '50x50':
            psffile = 'psf_mcao_msim_50mas.fits.fz' 

        #logging.info("Reading PSF: " + config_data["PSD_file"][AO_mode])
        ##psf_inp, head = fits.getdata(input_parameters["user_defined_psf"], 0, header=True, memmap=True)
        ##psf_wave = fits.getdata(input_parameters["user_defined_psf"], ext=1) 
        #psf_inp, head = fits.getdata(os.path.join(psf_path, config_data["PSD_file"][AO_mode]), 1, header=True, memmap=True)
        #psf_wave = fits.getdata(os.path.join(psf_path, config_data["PSD_file"][AO_mode]), 2)
        
        #pull the necessary file directly
        logging.info("Reading PSF: " + psffile)
        psf_inp, head = fits.getdata(os.path.join(psf_path, psffile), 1, header=True, memmap=True)
        psf_wave = fits.getdata(os.path.join(psf_path, psffile), 2)
        
        if psf_inp.ndim != 3:
            raise MSIMError("User PSF must have 3 dimensions.")
            
        if head['CDELT1'] != psfscale or head['CDELT2'] != psfscale:
            raise MSIMError("The PSF pixel scale must be = " + str(psfscale) + " mas.")
       
        user_psf = PchipInterpolator(psf_wave, psf_inp, extrapolate=True, axis=0)

    
    return user_psf


def create_psf(lamb, Airy=False):
    '''
    Returns a cube with the PSF for the given lambs generated from the PSD
 
    Args:
        lamb: lambda  [um]
        Airy: calculate Airy pattern
    
    Returns:
        cube: PSF
    '''
       
    global pup, stats, psd, xgrid_out, ygrid_out, jitter, psfscale, fov, diameter, AO_mode, rotation
    global zenith_seeing, air_mass
    print(AO_mode)
#    if AO_mode in ["MCAO", "AIRY"]:
#        # size of a pixel returned by psd_to_psf
#        psf_sampling = 2.
#        pix_psf = lamb*1e-6/(psf_sampling*diameter)*1/(4.85*1e-9) # mas
#        
#        if not Airy:
#            psf = psd_to_psf(psd, pup, diameter, phase_static = stats, lamb=lamb*1e-6, samp=psf_sampling, jitter=jitter/pix_psf)
#        else:
#            psf = psd_to_psf(psd*0., pup, diameter, phase_static = None, lamb=lamb*1e-6, samp=psf_sampling, jitter=np.repeat(0., 2))
#        
#        area_scale = (pix_psf/psfscale)**2
#        #print(area_scale)
#        if area_scale > 1:
#            # interpolate PSF
#            xgrid_in = (np.linspace(0, psf.shape[0]-1, psf.shape[0]) - psf.shape[0]*0.5)*pix_psf
#            ygrid_in = (np.linspace(0, psf.shape[1]-1, psf.shape[1]) - psf.shape[1]*0.5)*pix_psf
#            image = interp2d(xgrid_in, ygrid_in, psf, kind='cubic', fill_value=0.)
#            finalpsf = image(xgrid_out, ygrid_out)/area_scale
#        else:
#            # rebin PSF
#            side = int(psf.shape[0]*pix_psf/psfscale/2)*2
#            rebin_psf = frebin2d(psf, (side, side))/area_scale
#            center = side//2
#            finalpsf = rebin_psf[center-fov//2:center+fov//2, center-fov//2:center+fov//2]
#            
#        finalpsf[finalpsf < 0] = 0.
#        #fits.writeto("psf_orig.fits", psf, overwrite=True)
#        #print np.sum(finalpsf)
#        return finalpsf
# 
    if AO_mode in ["AIRY"]:
        # size of a pixel returned by psd_to_psf
        psf_sampling = 2.
        pix_psf = lamb*1e-6/(psf_sampling*diameter)*1/(4.85*1e-9) # mas
        
        psf = psd_to_psf(psd*0., pup, diameter, phase_static = None, lamb=lamb*1e-6, samp=psf_sampling, jitter=np.repeat(0., 2))
        
        area_scale = (pix_psf/psfscale)**2
        #print(area_scale)
        if area_scale > 1:
            # interpolate PSF
            xgrid_in = (np.linspace(0, psf.shape[0]-1, psf.shape[0]) - psf.shape[0]*0.5)*pix_psf
            ygrid_in = (np.linspace(0, psf.shape[1]-1, psf.shape[1]) - psf.shape[1]*0.5)*pix_psf
            image = interp2d(xgrid_in, ygrid_in, psf, kind='cubic', fill_value=0.)
            finalpsf = image(xgrid_out, ygrid_out)/area_scale
        else:
            # rebin PSF
            side = int(psf.shape[0]*pix_psf/psfscale/2)*2
            rebin_psf = frebin2d(psf, (side, side))/area_scale
            center = side//2
            finalpsf = rebin_psf[center-fov//2:center+fov//2, center-fov//2:center+fov//2]
            
        finalpsf[finalpsf < 0] = 0.
        return finalpsf
    
    elif AO_mode == "NOAO": # noAO Gaussian PSF
        # Beckers 1993 ARAA
        zenit_angle = np.arccos(1./air_mass)
        seeing_lambda = zenith_seeing/((lamb/0.5)**(1./5)*np.cos(zenit_angle)**(3./5))*1000. # mas
        sigma = seeing_lambda/2.35482
        
        Gauss2D = lambda x, y: 1.*np.exp(-(x**2 + y**2)/(2.*sigma**2))
        
        xx, yy = np.meshgrid(xgrid_out, ygrid_out)
        finalpsf = Gauss2D(xx, yy)
        finalpsf = finalpsf/np.sum(finalpsf)
        #fits.writeto("psf.fits", Gauss2D(xx, yy), overwrite=True)
        return finalpsf
    elif AO_mode == "MCAO":
        # user defined PSF
        return user_psf(lamb)
    #elif AO_mode == "MCI_LTAO" or AO_mode == "MCI_SCAO":
    #    return user_psf
    else:
        assert(1 == 0) # this should never happen

