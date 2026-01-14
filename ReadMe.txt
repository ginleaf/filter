+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
+ The Filament Trait-Evaluated Reconstruction (FilTER) code      +
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

This code is described in Panopoulou, Tassis, Goldsmith & Heyer, 2014, MNRAS,
444 (1): 2507-2524 (arXiv  1408.4809). 

It is distributed under the BSD-2 clause license, which is meant to
let you feel free to modify as you wish and send comments. 
Please cite the above in any related publication.

The code is still under development.
=======================
Requirements: 
(If you have troubles with other versions please contact panopg@physics.uoc.gr)
* python 2.7
* numpy  1.9.0
* pyfits 3.1.2
* scipy 0.9.0
* matplotlib 1.1.1rc
* pywcsgrid2
* DisPerSe (to provide a skeleton)
=======================

=======================
To run simply append the path to FilTER to your $PYTHONPATH.
Alternatively, add these lines to the script you are running:
import sys
sys.path.append(PATH_TO_FILTER)
Adjust parameters in main.py according to guidelines below.
Then run main.py .

Before running check that you can:
>>> import FilTER
>>> import pywcsgrid2 (and all above packages)
Have a DisPerSe output file of type .a.NDskl

To test the code (run as is), download this file from the Herschel Gould Belt Survey
archives (Andre et al. 2010), and save it in the sample folder:
www.herschel.fr/cea/gouldbelt/en/archives/taurusN3/taurusN3-500.fits.gz
=======================

=======================
Code description:

1) main.py: Parameters are defined in this script. Runs all necessary functions to 
         produce a set of filaments that are well-defined and have measured
         properties (e.g. width). Stores these in a numpy array. 
         Also calls some plotting functions.

Tips: Before running on the whole Disperse skeleton you can select to run the 
filstats.profiling routine (line 57 of main.py) for a small number of filaments 
e.g.  filstats.profiling(fils[:2], data, imsize, prof_ext, min_ext, noise_level, \
                       peakoff, plot_all, pixel2pc, units, threshold, binstep,\
                       pskip,fop, savedir)
This will help choose the parameters needed according to your own data.
Suggestions for parameters:

    - prof_ext = size of HALF profile in pixels i.e. how far away from the 
    spine of the filament you want the radial profile to reach. A good guess for this
    can be chosen by measuring by eye the radial size of filaments in your data (in pixels)
    and taking several times that value.

    - min_ext = dynamical fitting stops when profile reaches this size (pixels): 
    Choose this so that there are enough pixels to sample a significant part of the 
    peak of the profile
 
    - noise_level = threshold of intensity (whatever unit your image has) below which 
    peaks are discarded
    The noise level is easy to spot on most raw images.

    - plot_all = False # if True plots of all profiles will be made and saved. 
    Seriously slows down running time. Recommended only for testing if parameters 
    produce correct result in fitting.

    - distance = distance to cloud (needed for converting values in pixels to values in parsecc)

    - beamsize = size of beam in arcsec, for deconvolution of widths.
 
    - units = 'pc' No other type of units is supported currently.

    The following parameters have been optimized for a number of different data, so 
    I expect them to work in most data sets. Nevertheless, if something seems bad 
    (e.g. problems with the width estimation in many profiles) try playing with these too.

    - threshold = length to width ratio for filament to be acceptable 
    Just a definition (see paper).

    - binstep = make fwhm bins that are separated by this number of pixels (needed for mode estimation)
    The choice of this parameter will affect the precision of the widths that you get. 
    Convert this number to parsecs and you will get a feeling for a lower bound of the precision.
    If this is too large then the mode will not be correctly computed because you are making bins that are
    too big, if it is too small there will be problems with finding a mode (not enough bins) and consequently,
    defining the width of a profile.

    - peakoff = how far away can the peak of the gaussian fit be from the profile center (in pixels)?
    Depends on how strict you want to be with your selection criteria. This should be selected so that 
    the profile is centered on the correct filament and not a neighbouring one.

    - pskip = how many consecutive profiles can be 'bad'? Again, depends on how strict
    you want to be with your definition of 'continuity' of a structure.

2) filehandling.py: Includes functions to handle a.NDskl files (Disperse output), 
   al well as .fits files. Also includes Filament class.

3) filstats.py: Cuts profiles, fits them and assess them. Then new filaments are
   reconstructed according to the 3 criteria described in the paper.

4) prettyplots.py: Contains plotting functions to make the various figures in the paper.

5) sample: A sample DisPerSe skeleton of the taurus N3 image from the Herschel Gould Belt archives [1] 
is provided for testing the code. Parameters in main.py are given for this sample case.

6) run_disperse.py: An example of running disperse to create a skeleton and plot it.
=======================

Do not hesitate to ask questions or propose improvements!

Contact: georgia.panopoulou@chalmers.se

[1] http://www.herschel.fr/cea/gouldbelt/en/Phocea/Vie_des_labos/Ast/ast_visu.php?id_ast=66
