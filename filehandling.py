import numpy as np
import pyfits
import math

"""
Includes functions to handle a.NDskl files (Disperse output), 
al well as .fits files. Also includes Filament class.
"""
def read_fits(myfitsfile):
    """ 
    Parameters:
    (string) Path to 2-d fits file 

    Returns:
        (list) the size of the axes in pixels
    """
    hdulist=pyfits.open(myfitsfile)
    prihdr = hdulist[0].header
    imsize = [pyfits.getval(myfitsfile,'naxis1'),pyfits.getval(myfitsfile,'naxis2')]
    return imsize

def makenames():
    """
    Returns:
        list of useful names of fields in NDskl_ascii file
    """
    namescp2d = ['critical index','position0','position1','value','pairid','bound']
    namescp3d = ['critical index','position0','position1','position2','value','pairid','bound']
    namesfil2d = ['index cp1','index cp2','nsampl']
    namesfil3d = namesfil2d
    return [namescp2d,namescp3d,namesfil2d,namesfil3d]

def get_point_info(f,names,data = False):
    """
    Parameters:
    f: NDskl-ascii file handler
    names: list of strings with the type of data available for each critical point
           e.g. names = ['critical index','position0','position1','value','pairid','bound']
    data: (bool) set to True if you want to read CP_DATA section of file.
    
    Returns:
    (dictionary) with data in the CRITICAL POINTS section of the .a.NDskl file. 
    """
    cpoint = {}
    cols = f.readline().strip().split()
    assert len(cols) == len(names)
    for ii in range(0,len(cols)):
        cpoint[names[ii]] = float(cols[ii])
    if data == True:
        return cpoint
    cpoint['nfil'] = int(f.readline().strip())
    cpfil = []
    for ff in range(0,cpoint['nfil']):
        cpfil.append(map(int,f.readline().strip().split())) #converts list items into integers
    cpoint['cpfils'] = cpfil
    return cpoint

def get_fil_info(f,names):
    """
    Parameters:
    f: NDskl-ascii file handler
    names: list of strings with the type of data available for each filament, e.g. ['index cp1','index cp2','nsampl']

    Returns:
    Dictionary with data in the FILAMENTS section of the .a.NDskl file. 
    """
    fil = {}
    cols = f.readline().strip().split()
    assert len(cols) == len(names)
    for ii in range(0,len(cols)):
        fil[names[ii]] = int(cols[ii])
    sample_point_position = []
    for sample_point in range(0,fil[names[len(cols)-1]]):
        sample_point_position.append(map(float,f.readline().strip().split()))
    fil['sample point position'] = sample_point_position
    return fil

def get_fil_data_info(f,names,nsampl):
    """
    Parameters:
    
    f: (file handler) 
    names: (list) names of fields under FILAMENT section in NDskl-ascii file
    nsampl: (int) number of sampling points of current filament

    Returns:
    (list) List of sample point data of filament in the FILAMENTS DATA section of the a.NDskl file. 
           Form: [{filament's 1st sample point info},...{filament's nth sample point info}]."""
    fil = []
    sp = {}
    for sample_point in range(0,nsampl):
        cols = f.readline().strip().split()
        for ii in range(0,len(cols)):
            sp[names[ii]] = float(cols[ii])
        fil.append(sp)
    return fil

def read_NDskl_ascii(filename,namescp,namesfil):
    """
    Reads .a.NDskl file and returns dict with data for critical points and filaments.

    Parameters:

    filename: (string) NDSkl-ascii file as output by DisPerSe. Contains all info on image skeleton.
    namescp: (list of strings) Names of fields regarding Critical Points in NDskl_ascii file as returned by makenames()
    namesfil: (list of strings) Names of fields regarding Filaments in NDskl_ascii file as returned by makenames()

    Returns:

    (dictionary) Contains information on Critical Points and Filaments in a more code-friendly format than the NDskl-ascii file.

    Example: (for 2d case)
    >>>names = makenames()
    >>>section = read_NDskl_ascii(filename,names[0],names[2])
    """
    f = open(filename,'r')
    cpoints = []
    section = {}
    section['header'] = f.readline().strip()
    section['ndims'] = f.readline().strip()
    section['comments'] = f.readline().strip()
    section['bbox'] = f.readline().strip()
    section[f.readline().strip()] = cpoints
    cpoints.append(int(f.readline().strip()))
    for point in range(0,cpoints[0]):
        cpoints.append(get_point_info(f,namescp))

    fils = []
    section[f.readline().strip()] = fils
    fils.append(int(f.readline().strip()))
    for fil in range(0,fils[0]):
        fils.append(get_fil_info(f,namesfil))

    cpoints_data = []
    section[f.readline().strip()] = cpoints_data
    nfields = int(f.readline().strip())
    fields = []
    for field in range(0,nfields):
        fields.append(f.readline().strip())
    for cp in range(0,cpoints[0]):
            cpoints_data.append(get_point_info(f,fields,data = True))

    fils_data = []
    section[f.readline().strip()] = fils_data
    nfields = int(f.readline().strip())
    fields = []
    for field in range(0,nfields):
        fields.append(f.readline().strip())
    for filament in fils[1:]:
        fils_data.append(get_fil_data_info(f,fields,filament['nsampl']))

    f.close()
    return section

def angle_wrt_y(A,B):
    """
    Returns the angle between the line connecting two points in 2D and the positive y-axis.
    
    Parameters:
    A: (list) First point [x1,y1].
    B: (list) Second point [x2,y2].
    However, numpy default thinks [y,x] so keep that in mind.

    Returns:
    (float) Angle in range -pi/2, pi/2.
    """
    ax, ay = A
    bx, by = B
    if (by-ay) != 0.:
        return (180/math.pi)*math.atan((bx-ax)/(by-ay))
    else:
        return 90.0

class Filament(object):
    """
    The Filament class holds all information on filament objects and methods to analyse this information.
    """
    
    def __init__(self, points, image_size):
        """
        Parameters:

        points: (list) Contains pixels in the image that make up the Filament. Format is [[x1,y1],...,[xn,yn]].
        image_size: (list) The shape of the image as returned by numpy.
        """
        self.points = np.array(points)
        self.fwhms = np.zeros(len(points))
        
        assert type(image_size) == list
        assert len(image_size) == 2
        self.image_size = image_size

    def bounding_box(self):
        """
        Calculate the minimum bounding box.
        Return an array where first row is the top_left corner
        and second row is the bottom right corner of the box.
        """
        top_left = [self.points[:,0].min(), self.points[:,1].min()]
        bottom_right = [self.points[:,0].max(), self.points[:,1].max()]
        
        return np.array([top_left, bottom_right])
        
    def length(self):
        """
        Returns the length of the line that goes through all sampling points
        """
        last_point = self.points[0]
        length = 0
        for next_point in self.points[1:]:
            length = length + np.linalg.norm(next_point - last_point)
            last_point = next_point
        return length

    def angle(self):
        angles = []
        pointlist = self.points
        for pp in range(0,len(pointlist)-1): #check explicitly that closest point in pointlist is the one next to it
            angles.append(angle_wrt_y(pointlist[pp], pointlist[pp+1]))#equivalent that doesn't work: find_closest(pointlist[pp], pointlist[pp:])))
        return angles
    
    def median_width(self):
        # mask out zeros and nans before calculating the mean
        mask = np.logical_or(self.fwhms == 0., self.fwhms != self.fwhms)
        cleanfwhm = np.ma.array(self.fwhms,mask = mask)
        # compute median by accesing only valid entries of masked array
        return np.median(cleanfwhm[~cleanfwhm.mask])

    def mean_width(self):
        # mask out zeros and nans before calculating the mean
        mask = np.logical_or(self.fwhms == 0., self.fwhms != self.fwhms)
        cleanfwhm = np.ma.array(self.fwhms,mask = mask)
        # compute mean by accesing only valid entries of masked array
        return np.mean(cleanfwhm[~cleanfwhm.mask])

    def __repr__(self):
        return "Filament({0} points)".format(self.points.shape[0])

    def __str__(self):
        return self.__repr__()


def extract_filaments(ndskl_ascii_weird_format, image_size):

    filament_points = [ np.array(record['sample point position']) for record in ndskl_ascii_weird_format]  
    return [ Filament(points, image_size) for points in filament_points]

def load_filaments(filename, image_size):
    names = makenames();
    return extract_filaments(read_NDskl_ascii(filename, names[0], names[2])['[FILAMENTS]'][1:], image_size)

def fil_segs(filament,image_size): 
    # Filament segments as defined from sampling points become filament objects
    segs = [] 
    # deprecated: First point has to be treated differently (if it were a vector it would point from sampling point 1 to sampling point 0)
    # segs.append(Filament(np.array([filament.points[1],filament.points[0]]),image_size))
    # The last sampling point will be ignored
    for index in range(1,len(filament.points)):
        segs.append(Filament(np.array([filament.points[index-1],filament.points[index]]),image_size))
    return segs

