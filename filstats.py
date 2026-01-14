import math
import os
import numpy as np
import filehandling as fh
from scipy import optimize
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams.update({'font.size': 20, 'font.family': 'sans'})
'''
Module to cut across filament and produce intensity profile.
Has been developed to handle DisPerSe output.
'''
def Log(string, f):
    """ Write log messages
    """
    f.write(string)
    return

def quadratic(a,b,c):
        """solves ax**2+bx+c=0"""
        delta = b*b - 4*a*c
        assert delta >= 0.
        x1 = (-b + math.sqrt(delta))/(2*a)
        x2 = (-b - math.sqrt(delta))/(2*a)
        return x1, x2

def find_center(segment):
        """Takes Filament object as input.
           Finds center of segment and returns 
           pixel coordinates
        """
        x1, y1 = segment.points[0][0], segment.points[0][1]
        x2, y2 = segment.points[1][0], segment.points[1][1]
        x = abs(x1 - x2)/2. + min(x1,x2)
        y = abs(y1 - y2)/2. + min(y1,y2)
        return [x,y]

class Vect:
    #https://stackoverflow.com/questions/53970131/how-to-find-the-clockwise-angle-between-two-vectors-in-python
    def __init__(self, a, b):
        # initialize vector as Vect(math.cos(angle),math.sin(angle)),angle = angle of vector wrt x axis
        self.a = a
        self.b = b

    def findClockwiseAngle(self, other):
        # using cross-product formula
        inasin = (self.a * other.b - self.b * other.a)/(self.length()*other.length())
        # Catch a math domain error (e.g. if -1.000000001 in asin, set to -1)
        if  inasin < -1.:
            inasin = -1.
        if inasin > 1.:
            inasin = .1
        return -math.degrees(math.asin(inasin))
        
    def length(self):
        return math.sqrt(self.a**2 + self.b**2)


def perp_seg(segment,length):
        """
        Given a line segment of class Filament and a length
        returns the two end points of a line segment passing 
        through the first one's end point and whose
        bisector is the given segment
        """
        A,B = segment.points[0],segment.points[1]
        vect = [B[0]-A[0],B[1]-A[1]]
        mag = math.sqrt(vect[0]*vect[0]+vect[1]*vect[1])
        vect[0] = vect[0]/mag
        vect[1] = vect[1]/mag
        temp = vect[0]
        vect[0] = vect[1]
        vect[1] = -temp
        C = [B[0]+vect[0]*length,B[1]+vect[1]*length]
        D = [B[0]+vect[0]*(-length),B[1]+vect[1]*(-length)]
        return C,D

def perp_seg_end_always_on_same_side_of_fil(segment,perpsegment):
    # segment is along filament, perpsegment is perpendicular to it
    # Segment given must always have same direction (towards end of filament)
    A,B = segment.points[0],segment.points[1]
    # init Vect object for this segment
    # the cosine of the angle wrt x axis is x/length of segment
    input_segment_length = np.sqrt( (B[0]-A[0])**2 + (B[1]-A[1])**2)
    cosine = B[0] - A[0]
    sine = B[1] - A[1]
    filament_vect = Vect(cosine,sine)
    # Now do the same for the perp segment. Only for this one we don't know the 'right' direction
    C,D = perpsegment.points[0],perpsegment.points[1]
    # init Vect object for this segment
    # the cosine of the angle wrt x axis is x/length of segment
    input_segment_length = np.sqrt( (C[0]-D[0])**2 + (C[1]-D[1])**2)
    cosine = C[0] - D[0]
    sine = C[1] - D[1]
    perp_vect1 = Vect(cosine,sine)
    cosine = D[0] - C[0]
    sine = D[1] - C[1]
    perp_vect2 = Vect(cosine,sine)
    # Find angle that filament_vec makes with both of these vectors
    theta1 = filament_vect.findClockwiseAngle(perp_vect1)
    theta2 = filament_vect.findClockwiseAngle(perp_vect2)
    # Assign the vector with the positive angle to always define the start point of the segment
    if theta1 >= 0:
        return D,C
    elif theta2 >= 0:
        return C,D
    else:
        raise ValueError('Check angles of segments!')

def gaussian(B,x):
    ''' Returns the gaussian function for B = [mean,stdev,max] (no offset)
    '''
    return B[2]*np.exp(-((x-B[0])**2/(2*B[1]**2)))
 
def errfunc(B, x, y):
    return y-gaussian(B,x)

def fit_gauss(r_arr, y_arr, Ec, offset):
    ''' Fit a gaussian to the data.
        r_arr = distance values (x-axis)
        y_arr = data (y-axis)
        Ec = initial guess for gaussian peak value
        offset = background value. This will be subrtacted from the data before fitting
        Returns results of fit: peak position (mean), FWHM, values of the gaussian\
                with the fitted parameters, stdev, peak value (including offset)
    '''
    normdata = y_arr - offset
    p0 = [0., 1., Ec-offset]
    fit, covar, info, errmsg, ier = optimize.leastsq(errfunc, p0, \
                                    args=(r_arr, normdata), full_output = 1)
    #http://comments.gmane.org/gmane.comp.python.scientific.user/33879
    # Generate values from a gaussian with the fitted parameters (add offset)
    y2 = gaussian(fit, r_arr) + offset
    # The fit result for the mean of the gaussian is:
    gpeak = fit[0] 
    # Compute FWHM
    fwhm = 2.*math.sqrt(2.*math.log(2))*fit[1]
    #reduced_chi_square = ((y2-offset)**2).sum() / (len(y_arr) - len(p0))
    return gpeak, fwhm, y2, fit[1], fit[2]+offset

def sortit(xarr,yarr):
    ''' Input: xarr - array of values you want sorted
               yarr - array of values you want sorted according to sorted xarr
        Return: sorted xarr, yarr in order of xarr sorted
    '''
    return np.sort(xarr), yarr[np.argsort(xarr)]

def prof_statistics(gpeaks,fwhms,binned, max_ord=1):
    ''' 
        For a given profile calculate the mode FWHM, its frequency and how far the Gaussian fit is from the center 
        of the profile.

        Input:
            gpeaks = list of the position of the peak of all gaussian fits to profile (sorted according to fwhms)
            fwhms = list of fwhm of gaussian fits
            binned = output of numpy.histogram(fwhms). First item is array of frequencies, second is array of bin edges.
            max_ord = integer. 1 if you are checking for the first mode, 2 if checking for the second (second most frequent value in binned array)
        Returns:
            frequency of mode FWHM (float)
            position of gaussian fit corresponding to mode i.e. how far away it is from the center of the profile (int)
            FWHM (float)
    '''
    # Check if frequency array only has one item (only one bin of FWHM for this profile).
    if len(binned[0]) > 1:
        # sort frequencies (lowest to highest)
        sorted_indices = np.argsort(binned[0])
        # find highest frequency (or next to highest, depending on max_ord value)
        mode = binned[0][sorted_indices[-max_ord]]
        # find where the highest frequency is in the array
        #iimode = list(binned[0]).index(mode)
        iimode = sorted_indices[-max_ord]
        # count how many gaussian fits up to highest frequency
        nfitsmod = binned[0][:iimode].sum()
        # find the position of the peak of the gaussian that corresponds to the mode
        condition = abs(gpeaks[nfitsmod])
        # find the fwhm of the gaussian that corresponds to the mode
        fwhm = fwhms[nfitsmod]
    else:
        # if only one bin, the middle value is selected to be the FWHM
        # the middle value for an odd number of fwhms is different than for even (separating the two cases has no obvious effect)
        if len(fwhms)%2: # if odd select middle
            mid = len(fwhms)/2
        else: # if even, middle is between items, so select item on the lower side of middle
            mid = len(fwhms)/2 -1 
        mode = fwhms[mid]
        # the position of the peak of the middle gaussian will be checked
        condition = abs(gpeaks[mid])
        # fwhm assigned is that of the gaussian in the middle
        fwhm = mode
    return mode, condition, fwhm

def is_good_profile(gpeaks,fwhms,binned,prof_vals,noiselevel,peakoff):
    """ Check if profile is good: gaussian peak within peakoff pixels of
        the center, fwhm distribution has peak and peak above the
        noise level (set in the parameters section). 

        Input:
            gpeaks = list of the position of the peak of all gaussian fits to profile (sorted according to fwhms)
            fwhms = list of fwhm of gaussian fits
            binned = output of numpy.histogram(fwhms). First item is array of frequencies, second is array of bin edges.
            prof_vals = list of intensities along profile (nearest to center first, farthest last)
            noiselevel = level of noise intensity defined in main
            peakoff = maximum acceptable distance of the peak of a gaussian fit from the profile center (in pixels)
            
        Returns:
            True if profile is acceptable, False if not (bool)
            fwhm: the fwhm to be assigned to the profile (float)
    """
    mode, condition, fwhm = prof_statistics(gpeaks,fwhms,binned)
    # In case sth in the fit went wrong and you get a negative fwhm returned:
    if fwhm < 0.:
        fwhm = 0
        return False, fwhm
    # Check if profile has a mode FWHM, if it corresponds to a gaussian near the profile center and 
    # if the peak of the profile (max of points within peakoff of the center) is above the noise level
    if mode > 1 and condition < peakoff and np.max(prof_vals[:peakoff]) > noiselevel:
        return True, fwhm
    max_ord = 2 
    # Give a second chance, check for another mode
    mode, condition, fwhm = prof_statistics(gpeaks,fwhms,binned, max_ord = max_ord)
    if mode > 1 and condition < peakoff and np.max(prof_vals[:peakoff]) > noiselevel:
        return True, fwhm
    # If all else failed, this is a bad profile
    return False, 0

def extract_profiles(data, profiles, pixel2pc, units, binstep, extent, minextent,\
                    noiselevel, peakoff, counting, plot_all, beam, fop, psavedir):
    """Find filament width at each sampling point by dynamical fitting
       Input: 
       data = image (2d numpy array)
       profiles =  1d array of filaments profiles (Filament objects)
       pixel2pc = pixel to parsec conversion (float)
       units = units of distance (string)
       binstep = bin width in FWHM histogram (integer)
       extent = array of distances from 0 up to maximum extent of profile, i.e. prof_ext
       minextent = dynamical fitting ends when profile has reached this size in 
                   pixels (integer)
       noiselevel = level below which profiles are discarded
       peakoff = distance of gaussian fit peak from profile center in pixels 
       counting = filament id (integer)
       plot_all = if True plots of profiles will be made (boolean)
       beam = beam size in pc
       fop = logfile handle
       psavedir = folder to save profile images in

       Returns: 
       filament id, list of widths of profiles, list of profile assesement flags
    """
    posneg = np.array([])
    FWHM = np.array([])
    FWHM_deconv = np.array([])
    allfits = []
    fflags = np.array([])
    # Make the radial axis (on either side of the filament)
    for rr in extent:
        posneg = np.append(posneg, rr)
        # If you are on the spine, only append one point
        if rr != 0.:
            posneg = np.append(posneg, -rr)
    # Loop over filament profiles and get the intensity values along each profile
    for ii,profile in enumerate(profiles):
        prof_vals = np.array([])
        for pp,point in enumerate(profile.points):
            # Add pixel intensity value to the list of profile pixels
            try:
                prof_vals = np.append(prof_vals, data[int(point[1])][int(point[0])])
            except IndexError:
                Log('Profile outside image bounds - setting value = 0',fop)
                prof_vals = np.append(prof_vals, 0.)
        # Find the background of this profile. This will be subtracted from the
        # intensity values right before fitting.
        # Define a number of bins for intensities. Each bin will have at least 5
        # pixel values or 1/5 of the number of pixels of the profile
        intensity_bins = min(5, int(len(prof_vals)/5))
        # Make a histogram of intensities and select the background offset to be
        # the value of the minimum bin
        offset = np.histogram(prof_vals, intensity_bins)[1][0]
        # Dynamical profile. Will be changing.
        prof_valsd = prof_vals
        posnegd = posneg
        m = np.where(prof_vals == prof_vals.max())[0]
        Ec = prof_vals[0] 
        # Define lists to keep fit results from dynamical fitting of this profile
        fwhms = np.array([])
        gpeaks = np.array([])
        Cgaus = []
        ijk = 0
        profile.prof_vals = prof_vals
        # While profile is larger than min_ext, keep making it smaller
        # and fitting gaussians to it
        while len(posnegd) > minextent:      
            # fit gaussian to points in current profile
            gpeak, fwhm, fgauss, sigma, height = \
                                    fit_gauss(posnegd,prof_valsd, Ec, offset)
            # remember fwhm is in pixels
            # Keep track of how many fits have been made
            ijk +=1
            # Reduce profile by 2 pixels on each side
            prof_valsd = prof_valsd[:len(prof_valsd)-4]
            posnegd = posnegd[:len(posnegd)-4]   
            # Add this fit's fwhm and peak value to lists
            fwhms = np.append(fwhms, fwhm) 
            gpeaks = np.append(gpeaks, gpeak)
            # Make a gaussian function with these peak and fwhm values
            cgauss = gaussian([gpeak, sigma, height-offset], np.sort(posneg))+offset
            Cgaus.append(cgauss)
        # Sort the fwhm and peaks of the gaussians from dynamical fitting    
        fwhms, gpeaks = sortit(fwhms,gpeaks)
        # Make a distribution of the fwhm to find the mode value
        # If the profile passes quality checks, the mode found from this distribution
        # will be assigned as the width of this profile
        bins = int((fwhms.max()-fwhms.min())/binstep)
        # if the range of widths is larger than the binsize, then make bins with
        # the selected binsize (= binstep)
        if bins > 0: 
            binned = np.histogram(fwhms,bins = bins)
        # otherwise, number of bins is 1 (this sets the resolution of the algorithm to binstep)
        else:
            binned = np.histogram(fwhms, bins = 1)
        if plot_all: # Plot histograms of fwhm returned by dynamical fitting
            ff = plt.figure()
            axx = ff.add_subplot(111)
            if bins > 0:
                axx.hist(fwhms,bins = bins)
            else:
                axx.hist(fwhms)
            axx.set_ylabel('Number of fitted functions')
            axx.set_xlabel('FWHM of fitted functions (pixels)')
            plt.savefig(psavedir+'/'+'fil'+ str(counting) + 'prof'+ str(ii) +'fwhms.png')
            plt.clf()
            plt.close(ff) 

        Log('@ %d\n'%ii,fop)
        # Call is_good_profile to check whether profile fits our criteria (defined in main.py)
        is_good, fwhm = is_good_profile(gpeaks,fwhms,binned,prof_vals,\
                                            noiselevel,peakoff)  
        if is_good:
            FWHM = np.append(FWHM, fwhm)
            # Deconvolve the width with the beam size
            dfact = np.power(fwhm, 2) - np.power(beam/pixel2pc,2)
            if dfact > 0:
                FWHM_deconv = np.append(FWHM_deconv, np.sqrt(dfact))
            else:
                FWHM_deconv = np.append(FWHM_deconv, fwhm)
            # Assign a flag to this profile
            fflags = np.append(fflags, 0)
            if plot_all:
                f = plt.figure(2)
                axp = f.add_subplot(111)
                for cg in Cgaus:
                    axp.plot(np.sort(posneg)*pixel2pc,cg,'k--')
                axp.scatter(posneg*pixel2pc,prof_vals,color = 'y')
                axp.scatter(posnegd*pixel2pc,prof_valsd, color = 'b')
                sortr, sortf = sortit(posnegd*pixel2pc,fgauss)
                allfits.append(sortf)
                axp.plot(sortr,sortf,'c-') #sorted posneg and fgauss so plotted lines are correct
                axp.set_ylabel('Intensity')
                axp.set_xlabel('Distance from bone central axis (%s)'%units)
                plt.savefig(psavedir+'/'+'fil'+ str(counting) + 'prof'+ str(ii) +'.png')
                plt.clf()
                plt.close(2)

        else: 
            FWHM = np.append(FWHM, np.nan)  
            FWHM_deconv = np.append(FWHM_deconv, np.nan)
            # Assign a flag to the profile
            fflags = np.append(fflags,1)
        posnegd = posneg
        prof_valsd = prof_vals 
        fitting = True
    counting = counting + 1
    return counting, FWHM, FWHM_deconv, fflags


def convert_to_pc(filamentfile,resolution,savedir):
    filaments = np.load(savedir+filamentfile)
    for filament in filaments:
        filament.fwhms = list(np.array(filament.fwhms)*resolution)
        filament.length = list(np.array(filament.length)*resolution)
    filamentfile = filamentfile[:-4]+'parsecs.npy'
    np.save(savedir+filamentfile,filaments)
    return 

def newlengths(filaments, pskip, savedir, fop):
    """
    Divides filaments where there are more than pskip bad profiles
    Saves array of new filaments that only have good profiles (and less than pskip
    consecutive profiles).
    """
    newlendir = savedir+'/NewFilaments/'
    if not os.path.isdir(newlendir):
        os.mkdir(newlendir)

    nii = 0 # New filament counter
    newfils = filaments 
    
    widthsi = []
    fwhmsm = []
    imsize = filaments[0].image_size
    # Define new filaments
    for ii,fil in enumerate(filaments):
        Log('initial filament %d newfil %d\n'%(ii,nii),fop)

        tocheck = fil # Old filament to be investigated
        newfils = np.delete(newfils,nii) # Delete it from the newfils list
        allprofinds = np.arange(len(fil.profiles)) # Indices of all profiles (plus first one)
        goodprofinds = allprofinds[np.where(fil.fflags==0)[0]] # Keep only good profiles
        if len(goodprofinds)!=0: # If there are any good profiles
            
            newfilpoints = []
            newfilfwhms = []
            newfildeconvwidths = []
            newfilwidths = []
            npoints = 0
            newfilprofs = []
            
            ## Loop over profile indices
            for jj in range(0,len(goodprofinds)):
                
                ## Add first point to new filament if it is not already there
                if list(tocheck.points[goodprofinds[jj]]) not in newfilpoints:
                    
                    newfilpoints.append(list(tocheck.points[goodprofinds[jj]]))
                    npoints +=1
                    newfilfwhms.append(tocheck.fwhms[goodprofinds[jj]])
                    newfilprofs.append(tocheck.profiles[goodprofinds[jj]])
                    newfildeconvwidths.append(tocheck.deconvwidths[goodprofinds[jj]])
                
                try:
                   
                    ## Check if index of next good profile is within tolerance pskip
                    if goodprofinds[jj+1] <= goodprofinds[jj] + pskip:
                        newfilpoints.append(list(tocheck.points[goodprofinds[jj+1]]))
                        npoints+=1
                        newfilfwhms.append(tocheck.fwhms[goodprofinds[jj+1]])
                        newfilprofs.append(tocheck.profiles[goodprofinds[jj]])
                        newfildeconvwidths.append(tocheck.deconvwidths[goodprofinds[jj+1]])
                        
                    ## If not (and newfilpoints has more than 1 point), 
                    ## this is the end of the new filament. 
                    ## Add back to filament list.
                    else:
                        
                        if len(newfilpoints)>1: 
                            newfils = np.insert(newfils,nii,fh.Filament(newfilpoints,imsize))
                            newfils[nii].fwhms = newfilfwhms
                            newfils[nii].deconvwidths = np.array(newfildeconvwidths)
                            newfils[nii].profiles = newfilprofs
                            widthsi += newfildeconvwidths
                            fwhmsm.append(newfils[nii].median_width())
                            nii+=1  
                        newfilpoints = []
                        newfilfwhms =[]
                        newfilwidths = []
                        newfilprofs = []

                except IndexError: # Reached end of initial filament
                    
                    if len(newfilpoints) > 1:
                        Log('End fil, newpoints len %d \n'%len(newfilpoints),fop)
                        newfils = np.insert(newfils,nii,fh.Filament(newfilpoints,imsize))
                        newfils[nii].fwhms = np.array(newfilfwhms)
                        newfils[nii].deconvwidths = np.array(newfildeconvwidths)
                        newfils[nii].profiles = newfilprofs
                        widthsi += newfildeconvwidths
                        fwhmsm.append(newfils[nii].median_width())
                        nii+=1
                        
                    else: 
                        
                        Log('End fil, last point alone\n',fop)

    np.save(newlendir+'NewFilaments_pskip'+str(pskip)+'.npy',newfils)    
    np.save(newlendir+'FWHM.npy', np.array(widthsi))
    np.save(newlendir+'FWHMperfil.npy', np.array(fwhmsm))   
  
    return newfils

def profiling(fils, data, prof_ext, min_ext, noiselevel, peakoff,\
             pixel2pc, units, threshold, binstep, pskip, beam, fop, savedir, plot_all = False):
    """ Main function that calls all relevant functions to create profiles, 
        evaluate their properties (acceptable or not), find the width. 
        Saves into array of new filaments with only good profiles. 
        The 3:1 aspect ratio criterion is not included here. 
        ------
        Input: 
        fils = 1d numpy array of Filament objects 

        data = the image data as a 2d numpy array

        prof_ext = size of HALF profile in pixels (integer)

        min_ext = number of remaining pixels of profile to end dynamical fitting

        noiselevel = value of noise level (float)

        peakoff = max distance of peak of gaussian fit from profile center (pixels/integer)

        pixel2pc = size of the pixel in parsec

        units = distance units. If different than 'pc' then the pixel2pc 
                should be converted to your desired units (string)

        threshold = length to width ratio for filament to be acceptable (integer)

        binstep = make this many fwhm bins (needed for mode width estimation) (integer)

        pskip = if more consecutive 'bad' profiles then filament is split in two (integer)

        beam = size of telescope beam in pc. Will be used to deconvolve widths

        fop = logfile handler

        savedir = name of output directory

        plot_all = if True saves images of profiles in a separate folder (boolean)
        -------
        Output:
        saves arrays 'filaments_with_width.npy' containing Filament objects of the initial
        skeleton with widths assigned to their good profiles
        as well as 'NewFilaments_pskip*.npy' containing Filament objects that only contain
        acceptable profiles with widths assigned
    """
    imsize = list(data.shape[::-1])

    psavedir = savedir + '/profiles/'
    if not os.path.isdir(psavedir):
        os.mkdir(psavedir)

    extent = np.linspace(0, prof_ext, num=prof_ext+1)
    # Initialize array to save filaments and their widths 
    fils_with_width = []
    for counting, filament in enumerate(fils):
        Log('Currently at %d\n'%counting,fop)
        has_width = False
        profiles = [] # list of this filament's profiles
        psegs = [] # list of segments perpendicular to filament
        segments = fh.fil_segs(filament,imsize) # segments joining consecutive sampling points of a filament
        # Loop over all segments joining points along the filament spine
        for kk,seg in enumerate(segments):
            # Append 2 endpoints of segment perpendicular to seg 
            # choose length of 5 pixels, not important.. will only use direction
            perpendicular_segment_points = fh.Filament(perp_seg(seg,5), imsize)
            # Order the points so that the end point is always on the same side of the filament
            ordered_pseg_points = perp_seg_end_always_on_same_side_of_fil(seg,perpendicular_segment_points)
            pseg = fh.Filament(ordered_pseg_points, imsize)
            psegs.append(pseg)
            # Initialize points of filament profile at the position of this seg
            profile_pts = []
            # Find center point of perpendicular segment (center of profile)
            c = find_center(pseg)
            # Find dx and dy of endpoints of perpendicular segment
            dx = pseg.points[1][1]-pseg.points[0][1]
            dy = pseg.points[1][0]-pseg.points[0][0]
            # Loop over distances from the spine and keep the positions of points along the profile
            for r in extent:
                # Will find the two points that are at the intersection
                # of the line with slope of perpendicular segment and the circle 
                # with radius r, centered at the center of the profile (on the spine)

                # If perpendicular segment is parallel to x axis of image, do not let tangent go to inf
                if dy == 0:
                    # Find coordinates of the two points at distance r from the spine,  
                    # that lie on the line parallel to the x axis
                    x1, x2 = quadratic(1.,-2.*c[0],c[0]*c[0])
                    y1, y2 = c[1]+r,c[1]-r
                # For any other orientation of the segment:
                else:
                    # Calculate the slope of the perpendicular segment
                    tanf = dx/dy
                    # To solve circle equation (x-c0)**2 + (y-c1)**2 = r**2 (1)
                    # first substitute y-c1 from equation of line passing through [c0,c1], [x1,y1]
                    # i.e. y-c1 = tanf*(x-c[0]), to circle equation (1) and get:
                    # (x-c0)**2 (1+tanf**2) = r**2 (2)
                    # Calculate second parenthesis:
                    tat = 1. + tanf*tanf
                    # Solve (2) for x to get the x coordinates of the two points
                    x1, x2 = quadratic(1.,-2.*c[0],((c[0]*c[0])-r*r/tat))
                    # Get y of points from line equation: y = tanf*(x-c[0]) + c[1]
                    y1 = tanf*(x1-c[0]) + c[1]
                    y2 = tanf*(x2-c[0]) + c[1]

                # Only append one point to the profile if it is on the spine (r=0)
                if r == 0.:
                    profile_pts.append([x2,y2])
                else:
                    # To make sure the first point appended is on the same side of the filament
                    start_point,end_point = perp_seg_end_always_on_same_side_of_fil(seg,fh.Filament([[x1,y1],[x2,y2]], imsize))

                    profile_pts.append(end_point)#[x1,y1])
                    profile_pts.append(start_point)#[x2,y2])
            # Keep the points of the profile as a Filament object
            profiles.append(fh.Filament(profile_pts,imsize))
        # Make profiles, measure width and flag if good or bad
        counting, FWHMi, FWHM_deconv, fflags = extract_profiles(data, profiles, pixel2pc,\
                                  units, binstep, extent, min_ext, noiselevel, \
                                  peakoff, counting, plot_all, beam, fop, psavedir)
        
        # Add profiles to filament object
        filament.profiles = profiles
        # Add flags for profiles to filament object
        filament.fflags = fflags # profiles are flagged in this list (1 = bad)
        # There are filaments that no profile makes the test in extract_profiles
        # and FWHMi returns blanck
        if len(FWHMi) > 0:
            filament.fwhms = FWHMi
            # Deconvolved widths
            filament.deconvwidths = FWHM_deconv
            has_width = True        
            fils_with_width.append(filament)
    assert len(filament.fwhms) == len(filament.points)-1
    assert len(filament.fwhms) == len(filament.profiles)
    assert len(filament.fflags) == len(filament.profiles)
    np.save(savedir+'/filaments_with_width.npy',np.array(fils_with_width))

    newlengths(fils_with_width, pskip, savedir, fop)

    return

