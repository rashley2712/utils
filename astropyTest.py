#!/usr/bin/env python
import sys
try: 
	sys.path.remove('/home/astro/phsgan/Python64/lib/python/site-packages/astropy-0.3.2-py2.6-linux-x86_64.egg')
except (ValueError):
	print "No need to fix sys.path"
import astropy 
import photutils
import matplotlib.pyplot
import numpy as np
from photutils import datasets

print "Python version:", sys.version
print "Astropy version:", astropy.__version__

hdu = datasets.load_star_image()   
image = hdu.data[500:700, 500:700]   
image = hdu.data
print np.median(image)
image -= np.median(image)


from photutils import daofind
from astropy.stats import median_absolute_deviation as mad
bkg_sigma = 1.48 * mad(image)   
sources = daofind(image, fwhm=4.0, threshold=3*bkg_sigma)   
print sources

for s in sources:
	print s

figure = matplotlib.pyplot.figure(figsize=(10, 10))
matplotlib.pyplot.title("Sample image")
matplotlib.pyplot.imshow(image, cmap='gray')
#matplotlib.pyplot.gca().invert_yaxis()	
matplotlib.pyplot.show()

ax=figure.add_subplot(1,1,1)
matplotlib.pyplot.axis('off')

extent = ax.get_window_extent().transformed(figure.dpi_scale_trans.inverted())

figure.savefig('test.png', bbox_inches=extent)
