# Lensed Quasar Modeling
This repo will store the code I wrote to reconstruct images of gravitationally lensed quasars. This is important because it is difficult to extract good photometry from images. A good reconstruction will be one that is as good as an observation. In order to reconstruct an image of a lensed quasar there needs to be good estimates on the:
* positions of the quasar images
* position lensing galaxy (if not visible take average distance between the two quasars)
* flux, or amount of light given off by the quasars and lensing galaxy
* a measurement of the root mean square (RMS), and mean background signal of the images
* exposure time of the images

The code that will run create the lens models and reconstructions is ``LenstronomyMaster.py``. The code makes use of [lenstronomy](https://github.com/sibirrer/lenstronomy), a gravitational lensing software package.

# Requirements
* `CosmoHammer`
