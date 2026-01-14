import os
import pyfits
import filehandling as fh
import numpy as np
import filstats
import pretty_plots as pp

def log(string, f):
    """ Write log messages
    """
    f.write(string)
    return

if __name__ == '__main__':

    # Define file names and output folder
    # Your image/data file (.fits format only)
    fitsfile = 'sample/taurusN3-500.fits'
    # Pixel size (CDELT1) unit in the header
    cdeltu = 'degrees' # 'degrees' or 'arcmin' or 'arcsec'
    # The Disperse output in a.NDskl format (or equivalent numpy array '*.npy')
    skeletonfile = os.path.dirname(fitsfile)+'/DisperseSkeleton.a.NDskl'
    # The folder where the output of main will be stored
    savedir = os.path.dirname(skeletonfile)+'/savedir/'
    
    ### PARAMETERS ###
    prof_ext = 20 # size of HALF profile in pixels (Choose how far away from the 
                 # spine of the filament you want the radial profile to reach)
    min_ext = 6 # dynamical fitting stops when profile reaches this size (pixels)   
    noise_level = 10. # threshold of intensity (whatever unit your image has) below which peaks are discarded
    plot_all = False # if True plots of all profiles will be made and saved. Seriously slows down running time. Recommended only for testing if parameters produce correct result in fitting.
    distance = 140.# cloud distance in pc
    units = 'pc' # distance units
    beamsize = 36 # beam size in arcsec, for deconvolution of widths
    threshold = 3 # length to width ratio for filament to be acceptable 
    binstep = 5 # fwhm bin size in pixels (needed for mode width estimation)
    peakoff = 7 # how far away can the peak of the gaussian fit be from the profile center (in pixels)?
    pskip = 3 # how many consecutive profiles can be 'bad'?
    ##################
    
    # Read image and get size
    data = pyfits.getdata(fitsfile)
    imsize = list(data.shape[::-1])
    header = pyfits.getheader(fitsfile)
    # Calculate the pixel size in pc
    # pc to AU conversion
    pc2au = 206264.
    # arcsec conversion
    if cdeltu == 'degrees':
        asconv = 3600.
    if cdeltu == 'arcmin':
        asconv = 60.
    if cdeltu == 'arcsec':
        asconv = 1.
    pixel2pc = distance*np.absolute(header['CDELT1'])*asconv/pc2au
    # Beam size in pc
    beampc = distance*beamsize/pc2au
    # Get rid of nans as they may cause problems further on
    data = np.nan_to_num(data)
    # Check file type and load skeleton as an array of Filament objects
    if skeletonfile[-3:]=='skl':
        fils = fh.load_filaments(skeletonfile, imsize)
    else:
        if skeletonfile[-3:]=='npy':
            fils = np.load(skeletonfile)
        else: 
            raise Exception('Wrong type of input file')
    
    # Create savedir if needed
    if not os.path.isdir(savedir):
        os.mkdir(savedir)
    fop = open(savedir+'status.log','w',0) # will flush the buffer after every write
    log('Parameters chosen: \n',fop)
    log('prof_ext = %d\n'%prof_ext,fop)
    log('min_ext = %d\n'%min_ext,fop)
    log('noise_level = %.1f\n'%noise_level,fop)
    log('pixel2pc = %.3f\n'%pixel2pc,fop)
    log('threshold = %d\n'%threshold,fop)
    log('binstep = %d\n'%binstep,fop)
    log('peakoff = %d\n'%peakoff,fop)
    log('pskip = %d\n'%pskip,fop)
    log('Starting to make filament profiles\n',fop)
    print 'Number of filaments in skeleton:',len(fils)
    print 'Starting to cut and fit profiles'

    # Run the main analysis
    filstats.profiling(fils, data, prof_ext, min_ext, noise_level, \
                       peakoff, pixel2pc, units, threshold, binstep,\
                       pskip, beampc, fop, savedir, plot_all = plot_all)
    
    # Make plots
    # Load in all 'filaments' of the skeleton with width info
    allfils = np.load(savedir+'filaments_with_width.npy')
    # Make a distribution of widths of the entire skeleton (fig 9 right)
    pp.fil_FWHM_distro(allfils,pixel2pc,savedir) 
    # Same for widths that have been deconvolved with the beam size
    pp.fil_FWHM_distro(allfils,pixel2pc,savedir, deconv = True)
    # Load in filaments that comply to the quality checks (reprocessed skeleton)
    newfils = np.load(savedir+'NewFilaments/'+'NewFilaments_pskip'+str(pskip)+'.npy')
    # Make a list of filaments that have length to width ratio >= threshold
    true_fils = []
    filindices = []
    # Loop over filaments from the reprocessed skeleton and add the 'true' ones to the list
    for nn in range(len(newfils)):
        if newfils[nn].length() >= threshold * newfils[nn].mean_width(): #comparison is in pixels here
            true_fils.append(newfils[nn])
            filindices.append(nn)
    # Plot the 'true' filaments skeleton on the fits image
    pp.true_fil_plot(newfils, threshold, units, pixel2pc, savedir)
    # Plot distribution of filament widths
    # You can pass in true_fils instead of newfils to get the distribution of 
    # widths of only the those passing the ratio threshold
    pp.fil_FWHM_distro(newfils,pixel2pc,savedir) 
    # Loop over the 'true' filaments
    for trueindex, filindex in enumerate(filindices):
        # Plot mean profile of true filaments (fig 13)  
        pp.meanprof_linlog(true_fils[trueindex],filindex,pixel2pc,prof_ext,fitsfile,savedir)
    # Plot all filaments on image (fig 8)
    pp.plot_colored_fils(savedir,fils,fitsfile,savename='colored_fils.png')
    pp.plot_colored_fils(savedir+'/NewFilaments/',newfils,fitsfile,savename='colored_newfils.png')
    pp.plot_colored_fils(savedir+'/NewFilaments/',true_fils,fitsfile,savename='colored_truefils.png')

    log('The end. Output plots & files in %s\n'%savedir,fop)


