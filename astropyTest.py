#!/usr/bin/env python
import astropy 
import sys
import photutils
import numpy as np
from photutils import datasets

print "Python version:", sys.version
print "Astropy version:", astropy.__version__

hdu = datasets.load_star_image()   
image = hdu.data[500:700, 500:700]   
image -= np.median(image)

from photutils import daofind
from astropy.stats import median_absolute_deviation as mad
bkg_sigma = 1.48 * mad(image)   
sources = daofind(image, fwhm=4.0, threshold=3*bkg_sigma)   
print sources

