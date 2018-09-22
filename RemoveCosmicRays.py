"""
Given an image file (.fits format) this script will remove cosmic rays from image which will contaminate the data
The number of iterations the cleanup process will occur is user specified as wll as the various parameters of the clen up process such as gain, noise, and standard deviation of noise etc...
"""

import cosmics

"Import the image file"
file_path = "/Users/edenmolina/Documents/Astro/Coadds/2013/1206_10_1.fits"
name = "J1206"
image_file, hdr = cosmics.fromfits(file_path)
iterations = 10
"Create lacosmics object"
c = cosmics.cosmicsimage(image_file, gain=7, readnoise=.4, sigclip = 2, sigfrac = .5, objlim = 4)

"Clean the image and export"
for i in range(iterations):
    print "Iteration: %s" %i
    c.lacosmiciteration()
    c.clean()

cosmics.tofits('/Users/edenmolina/Desktop/Masks/%s_clean.fits'%name, c.cleanarray, hdr)
cosmics.tofits("/Users/edenmolina/Desktop/Masks/%s_mask.fits"%name, c.mask, hdr)