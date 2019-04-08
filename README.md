# Lensed Quasar Modeling
This repo will store the code I wrote to reconstruct images of gravitationally lensed quasars. This is important because it is difficult to extract good photometry from images. A good reconstruction will be one that is as good as an observation. In order to reconstruct an image of a lensed quasar there needs to be good estimates on the:
* positions of the quasar images
* position lensing galaxy (if not visible take average distance between the two quasars)
* flux, or amount of light given off by the quasars and lensing galaxy
* a measurement of the root mean square (RMS), and mean background signal of the images
* exposure time of the images

The code that will run create the lens models and reconstructions is `LenstronomyMaster.py`. The code makes use of [lenstronomy](https://github.com/sibirrer/lenstronomy), a gravitational lensing software package. The lenses are modeled as Single Isothermal Ellipsoids (SIE) with a shear component, the source of the quasars are modeled as Sersic profile, and the lensing galaxies-if visible- are modeled as Sersic ellipses. A mask can be input into the modeling in order to shroud stars or cosmic rays which will only confuse lenstrnonmy-as it might these these extra light sources are additional quasar images. The code will determine the correct combination of parameters to accurately describe the systems in the images. The solutions to the lensing models, and lensing equations are optimized using Particle Swarm Optimization (PSO) and Markov Chain Monte Carlo (MCMC).



# Requirements
* `Python 2.7`
* `CosmoHammer`
* `astropy`
* `numpy`
* `scipy`
* `matplotlib`
* `fastell4py`
