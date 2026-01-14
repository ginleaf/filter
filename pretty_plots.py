import numpy as np
import matplotlib.pyplot as plt
import math
import os,sys
import matplotlib
import pyfits
import scipy.stats
from matplotlib.font_manager import FontProperties
from matplotlib import rc, font_manager
import pywcsgrid2
import filehandling as fh
import warnings
warnings.simplefilter("ignore", Warning)

sizeOfFont = 8
fontProperties = {'family':'serif','serif':['Times New Roman'],
    'weight' : 'normal', 'size' : sizeOfFont}
ticks_font = font_manager.FontProperties(family='Times New Roman', style='normal',
    size=sizeOfFont, weight='normal', stretch='normal')
rc('text', usetex=True)
rc('font',**fontProperties)
rc('legend',frameon=False)

class supervised_plotting(object):

    def __init__(self, _total_colors=10000000):
        self.total_colors = _total_colors
        self.colors = np.random.rand(self.total_colors,3)
        self.color_index = 0
    
    def get_next_color_index(self):
        self.color_index += 1
        if self.color_index == self.total_colors:
            self.color_index = 0
        return self.color_index

    def get_next_color(self):
        rgb_color = self.colors[self.get_next_color_index()]
        rgb_color *= 255
        rgb_color = np.array(rgb_color, dtype="int")
        return "#{0:02X}{1:02X}{2:02X}".format(rgb_color[0], rgb_color[1], rgb_color[2])
    
    def plot_fil(self,fil, pl = True, invert_y = False,color = None,lw = 1.,axis = 0.):
        """plots sampling points of a Filament"""
        #same is achieved with zip(*fil['sample point position'])   
         
        if not axis:
            fig = plt.figure()
            ax = plt.gca()
        else:
            ax = axis
        x = []
        y = []
        for sp in fil.points:
                x.append(sp[0])
                if invert_y:
                    y.append(fil.image_size[1] - sp[1])
                else:
                    y.append(sp[1])
        if pl == True:
            if color == None:
                ax.plot(x,y,'-',lw =lw,c=self.get_next_color())            
            else:
                ax.plot(x,y,'-',lw =lw,c=color)                
                
        return [x,y]

def plot_colored_fils(savedir,fils,fitsfilename,
                      color = None,skycoords=False,savename='colored_truefils.png'):
    """ Makes plot like fig 8: overplot disperse filaments on data image"""
    fig = plt.figure(figsize=(7,4.5))
    ax = plt.gca()
    a=supervised_plotting()
    if not skycoords:
        data = np.flipud(pyfits.getdata(fitsfilename))
        dat = np.nan_to_num(data)
        plt.gray()
        plt.imshow(dat)
        plt.gca().xaxis.set_visible(False)
        plt.gca().yaxis.set_visible(False)
    else:
        f = pyfits.open(fitsfilename)
        h, dat = f[0].header, f[0].data
        ax = pywcsgrid2.subplot(111, header=h)  
        plt.gray()
        plt.imshow(1-dat,origin="lower")       
    plt.axis('image')
    fig.tight_layout()
    for filament in fils:
        if not skycoords:
            a.plot_fil(filament, invert_y = True, color=color,axis=ax) 
        else:
            a.plot_fil(filament, invert_y = False, color= color,axis=ax)
    plt.savefig(savedir+savename,bbox_inches='tight')
    plt.clf()
    plt.close(fig)
    return

def true_fil_plot(fils, threshold, units, pixel2pc, savedir):
    assert units == 'pc' # NOT PIXELS!
    
    true_fils = []
    lengths = []
    widths = []
    for kk,filament in enumerate(fils):
        length = filament.length() * pixel2pc
        width = filament.median_width()* pixel2pc
        ratio = length / width

        if ratio > threshold:
            true_fils.append(filament)
       
        lengths.append(length)
        units = units
        if width == 0:
            widths.append(np.inf)
        else:
            widths.append(width)
    widths = np.array(widths)
    lengths = np.array(lengths)
    true_fils = np.array(true_fils)
    np.save(savedir+'/true_filaments.npy',true_fils)
    condition = (lengths/widths) > threshold
    glengths = lengths[condition]
    gwidths = widths[condition]
    # Make plot like fig 9 (left)
    fig = plt.figure(1,figsize=(3.45,2.5))
    ax = plt.gca()
    linepoints = np.linspace(0, max(lengths), num = 3)
    line1 = linepoints
    line3 = linepoints/3
    ax.plot(linepoints,line1,'--',dashes = (3,4),color='0.5',label = 'Aspect ratio = 1')
    ax.plot(linepoints,line3,'-',color='0.5',label = 'Aspect ratio = 3')
    ax.scatter(lengths, widths, c='k', marker = 'o',s = 1,lw = 0.5)
    ax.scatter(glengths, gwidths, marker = 'o',c = 'r',s = 9, lw = 1)
    ax.set_xlabel('Filament lengths (pc)')
    ax.set_ylabel('Filament widths (pc)')
    left  = 0.15  # the left side of the subplots of the figure
    right = 0.9    # the right side of the subplots of the figure
    bottom = 0.13   # the bottom of the subplots of the figure
    top = 0.9      # the top of the subplots of the figure
    wspace = 0.4   # the amount of width reserved for blank space between subplots
    hspace = 0.3  # the amount of height reserved for white space between subplot    
    plt.subplots_adjust(left=left, bottom=bottom, right=None, 
                    top=None, wspace=wspace, hspace=hspace)
    plt.savefig(savedir+'/length_width.png')
    plt.close()
    return    

def fil_FWHM_distro(fils,pixel2pc,savedir, deconv = False, bins = 40):
    FWHMind = np.array([])
    for ii, filament in enumerate(fils):
        if deconv:
            # Add filament deconvolved widths to the list
            FWHMind = np.append(FWHMind, filament.deconvwidths)
            dstr = 'deconvolved'
        else:    
            # Add filament widths to the list (discarding nans)
            FWHMind = np.append(FWHMind,filament.fwhms)
            dstr = ''
    #mask any zeros or nans
    mask = np.logical_or(FWHMind == 0., FWHMind != FWHMind)
    cleanfwhm = np.ma.array(FWHMind, mask = mask)
    FWHMind = cleanfwhm[~cleanfwhm.mask] * pixel2pc
    stdev = np.std(FWHMind)
    print 'Stdev of distribution is',stdev, 'Mean is', np.mean(FWHMind)
    # Make the canvas
    fig = plt.figure(1)
    ax0 = fig.add_subplot(111)
    indbins = ax0.hist(FWHMind,bins = bins, histtype = 'step',color='k')
    ax0.set_xlabel('%s FWHM along filaments (pc)'%dstr)
    ax0.set_ylabel('Number of profiles')
    plt.savefig(savedir+'/' + 'FWHM-hist'+dstr+'.png')
    plt.clf()
    plt.close()
    return
    
def meanprof_linlog(newfilament,counting,pixel2pc,ext,fitsfile,psavedir,sigmasubtract=13,ecadd=1):
    """ Plots mean profiles like fig 13, 14"""
    # Radial extent of filament profiles (in pixels)
    extent = np.linspace(0, ext, num= ext+1)
    # Make the radial axis (on either side of the filament)
    posneg = [0.]
    for rr in extent[1:]:
        posneg.append(rr)
        posneg.append(-rr)
    posneg = np.array(posneg)
    # The data image size
    imsize = newfilament.image_size
    # Create canvas for plot
    f = plt.figure(1,figsize=(7,2.6))
    ax = f.add_subplot(121)  
    ax1 = f.add_subplot(122)
    # First make linear plot of mean profile
    data = pyfits.getdata(fitsfile)
    dat = np.nan_to_num(data)
    # Small inlet with filament
    b = plt.axes([0.34, 0.68, .15, .2])
    plt.imshow(1-dat)
    plt.gray()
    bb = list(newfilament.bounding_box())
    plt.axis([bb[0][0]-20,bb[1][0]+20,bb[0][1]-20,bb[1][1]+20])
    ssegs = fh.fil_segs(newfilament,imsize)
    for seg in ssegs:
        plt.plot(seg.points[0][0],seg.points[0][1],'w,',lw =0.1)
        plt.plot(seg.points[1][0],seg.points[1][1],'w,',lw =0.1)
    plt.setp(b,xticks=[],yticks=[])
    ax.axis([-extent[-1]*pixel2pc,extent[-1]*pixel2pc,-1e1,1e3])
    
    even = np.array([i for i in range(0,len(newfilament.profiles),4)])
    # This is because too many points make huge figure.. Optional if up to 200.
    fewer_profiles = np.array(newfilament.profiles)[even]
    # List that keeps all values of intensity of all profiles
    all_vals= []
    # Loop over (even) profiles of this filament and plot intensity values of linear profile
    for profile in fewer_profiles:
        ax.scatter(posneg*pixel2pc,profile.prof_vals,color='0.75',s=1)
        # For logarithmic plot you will only be plotting the positive radii of the profile
        profile.pos = posneg[posneg>0]
        pos_prof_vals = profile.prof_vals[posneg>0]
        all_vals.append(profile.prof_vals)
        ax1.scatter(np.log10(profile.pos*pixel2pc),np.log10(pos_prof_vals),color='0.75',s=1)
    # Fix axis properties for linear plot    
    ax.set_autoscaley_on(True)
    ax.set_xlim([-extent[-1]*pixel2pc,extent[-1]*pixel2pc])
    amin, amax = np.min(all_vals), np.max(all_vals)
    ax.set_ylim([amin,amax+amax/4])
    ax.set_ylabel('Integrated intensity (K km/s)')
    ax.set_xlabel('Distance from central axis (pc)')
    mean_vals = np.mean(all_vals,0)
    sortr, sortm = np.sort(posneg), mean_vals[np.argsort(posneg)]
    Ec = mean_vals[0]
    ax.plot(sortr*pixel2pc, sortm,'-k',linewidth = 2)
    # Small inlet with map
    a = plt.axes([0.07, 0.68, .15, .2])
    plt.imshow(1-dat)
    plt.gray()
    sp = supervised_plotting()
    sp.plot_fil(newfilament, color= 'k', axis=a)
    a.axis([0,imsize[0],0,imsize[1]])
    plt.setp(a,xticks=[],yticks=[])
    
    # Now set up logarithmic plot
    ax1.set_ylabel('log(Integrated intensity)')
    ax1.set_xlabel('log(Distance) from central axis ')
    ax1.plot(np.log10(sortr[sortr>0]*pixel2pc), np.log10(sortm[sortr>0]),'-k',linewidth = 1.5, label = 'Mean')
    ax1.legend(loc = 3,prop={'size':8})
    ax1.set_ylim(ymin=-3)
    ax1.set_xlim(xmax = np.log10(max(posneg*pixel2pc)))
    left  = 0.125  # the left side of the subplots of the figure
    right = 0.9    # the right side of the subplots of the figure
    bottom = 0.1   # the bottom of the subplots of the figure
    top = 0.9      # the top of the subplots of the figure
    wspace = 0.4   # the amount of width reserved for blank space between subplots
    hspace = 0.3  # the amount of height reserved for white space between subplot    
    plt.subplots_adjust(left=None, bottom=None, right=None, 
                    top=None, wspace=wspace, hspace=hspace)
    f.tight_layout()
    plt.savefig(psavedir+'fil'+ str(counting) + 'log-lin-light.png',bbox_inches='tight',dpi=300)
    plt.close(1)  
    return
