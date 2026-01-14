import subprocess
import os
import sys
import numpy as np
from pretty_plots import plot_colored_fils
from filehandling import load_filaments
import pyfits
# Run mse and then skelkonv on a fits image to get its skeleton
# Arguments: fits image, 0 or 1 (for using persistence cut value or choosing with pdview), persistence threshold value
# example: python sample/taurusN3-500.fits 0 3
fits, extension = os.path.splitext(sys.argv[1])
full_path_to_file = os.path.abspath(fits+extension)
fitsbasename = os.path.basename(fits)
workdir = os.getcwd()

print fitsbasename
persistence = int(sys.argv[2]) #0 for specifying persistence cut or 1 for showing pdview (interactive)
cut = sys.argv[3]

#If you want to check out the persistence diagram before making a threshold cut
if persistence:
    if os.path.isfile(fits+'.MSC'):
        subprocess.call(["mse", fits +".fits", "-periodicity","0", "-robustness","-interactive" ,"-upSkl", "-outName", fits, "-loadMSC", fits+".MSC"])
    else:
        subprocess.call(["mse", fits +".fits","-periodicity", "0","-robustness","-interactive","-upSkl", "-outName", fits])
    
    cut = raw_input('What threshold did you select (round to first decimal)?')
    folder = os.path.dirname(os.path.abspath(fits+extension))+'/'+'persist'+cut+'/'
    if not os.path.isdir(folder):
        os.mkdir(folder)
    os.rename(fits+'_c'+cut+'.up.NDskl', folder+fitsbasename+'_c'+cut+'.up.NDskl')
# Otherwise it will just cut where you told it to in the input
else:
    folder = os.path.dirname(os.path.abspath(fits+extension))+'/'+'persist'+cut+'/'
    if not os.path.isdir(folder):
        os.mkdir(folder)
    if os.path.isfile(fits+'.MSC'):
        subprocess.call(["mse", fits +".fits", "-periodicity","0", "-robustness","-cut", cut ,"-upSkl", "-outName", folder+fitsbasename , "-loadMSC", fits+".MSC"])
    else:
        subprocess.call(["mse", fits +".fits","-periodicity", "0","-robustness","-cut", cut,"-upSkl", "-outName", fits])
        os.rename(fits+'_c'+cut+'.up.NDskl', folder+fitsbasename+'_c'+cut+'.up.NDskl')

print 'Done running mse. Give me input for skelkonv please'
rob = raw_input('Robustness level is?')
asb = raw_input('Assemble at what angle (0-90)?') #say 70
smooth = raw_input('Smooth how much?') # say 30...

smoothstr = str(smooth)
if len(smoothstr) == 2:
    smothstr = '0'+smoothstr

#Now go to the resulting folder and run skelkonv there
os.chdir(folder)
print 'Went in ', folder
skel = fitsbasename+'_c'+cut+'.up.NDskl'
subprocess.call(["skelconv", skel, "-trimBelow", "0","-breakdown","-smooth", smooth,"-trimBelow", "robustness", rob, "-assemble", asb,"-to","NDskl_ascii","-outName", folder + 'c'+cut+'.up.NDskl' + ".TRIM0.BRK.S%s.TRIMrob"%smoothstr+rob+".ASMB"+asb,"-noTags"])
print 'Saved as ', folder + 'c'+cut+'.up.NDskl' + ".TRIM0.BRK.S%s.TRIMrob"%smoothstr+rob+".ASMB"+asb

# Now prepare for the plotting functions
os.chdir(workdir)
# Read image and get size
data = pyfits.getdata(full_path_to_file)
data = np.nan_to_num(data)
imsize = list(data.shape[::-1])
skeletonfile = folder + 'c'+cut+'.up.NDskl' + ".TRIM0.BRK.S%s.TRIMrob"%smoothstr+rob+".ASMB"+asb+'.a.NDskl'
fils = load_filaments(skeletonfile, imsize)

plot_colored_fils(folder,fils,full_path_to_file,
                      color = None,skycoords=False,savename='colored_truefils_rob'+rob+'.png')
