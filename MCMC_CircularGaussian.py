"""
Stable Version of the MCMC

From an image file the user will click on a part of the image and this scipt, using a Markov Chain Monte Carlo, will fit the flux (data values at the point) to a circular Gaussian
This code is a support function for the scipt named "grabCoordinates_v2.py".

"""

import matplotlib.pyplot as plt
import numpy as np
import astropy.io.fits as fits
import emcee
from scipy.stats import multivariate_normal
import scipy.optimize as op
import corner
from matplotlib.colors import LogNorm
import pandas as pd
from matplotlib.ticker import MaxNLocator
import multiprocessing
import sys
import os
import random

def SaveModel(data, Model, SavePath, flux):

    plt.imshow(data, cmap="gray", norm = LogNorm(), vmin = 1, vmax = 25, origin = 'lower')
    plt.contour(Model, origin='lower')
    plt.colorbar()
    #plt.show()
    plt.savefig(SavePath + "%sJan21_1950.png" % flux, dpi=50)
    plt.close()

"Liklihood function"
def lnlike(theta, x, y, data):
    flux, mu_x, mu_y, sigma_x = theta

    covariance = np.array([[sigma_x ** 2, 0], [0, sigma_x ** 2]])

    x_m, y_m = np.mgrid[0:len(x), 0:len(y)]
    pos = np.dstack((x_m, y_m))

    rv = multivariate_normal([mu_x, mu_y],covariance)
    Model = rv.pdf(pos)*flux

    SaveModelBool = False
    if SaveModelBool:
        SaveModel(data, Model, "/Users/edenmolina/Desktop/Test Export/", flux)
    return -0.5 * (np.sum((data - Model) ** 2))

"Prior probability"
def lnprior(theta):
    flux, mu_x, mu_y, sigma_x = theta

    bounds = 50
    if 1 < flux < 200000 and 10 < mu_x < 90 and 10 < mu_y < 90 and 0. < sigma_x < bounds:
        return 0.0
    else:
        return -np.inf

"Probaility Distribution"
def lnprob(theta, x, y, data):
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + lnlike(theta, x, y, data)

def SaveSamplerChainPlot(axes, axes_number, data, parameter_name):
    axes[axes_number].plot(data, color="k", alpha=0.2)
    axes[axes_number].yaxis.set_major_locator(MaxNLocator(5))
    axes[axes_number].set_ylabel(parameter_name)

def Flux(SampleName, iteration_number ,data, flux_c, mu_x_c, mu_y_c, mu_x_correction, mu_y_correction):
    "Load the Image Data"
    data = data

    # Number of Walkers and Rank of the parameter space
    nwalkers, ndim = 20, 4

    print "Image Peak", flux_c

    "Initial Guesses for the 6 parameters to Optimize"
    flux = np.random.normal(flux_c, 10000, nwalkers)
    for i in range(len(flux)):
        if flux[i] < 0:
            flux[i] = -flux[i]
    #flux = np.linspace(flux_c, 100000, nwalkers)
    mu_x = np.random.normal(mu_x_c, 10, nwalkers)
    mu_y = np.random.normal(mu_y_c, 10, nwalkers)
    sigma_x = np.random.normal(15, 10, nwalkers)
    for i in range(len(sigma_x)):
        if sigma_x[i] < 0:
            sigma_x[i] = -sigma_x[i]

    "Initializes points to plug in for position"
    x = np.linspace(0, len(data)+1, len(data))
    y = np.linspace(0, len(data[0])+1, len(data[0]))

    "Put Together the Initial Guesses for the Model"
    p0 = np.vstack((flux, mu_x, mu_y, sigma_x)).T

    # Shape the Gaussian list to an array to visualize the data
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(x, y, data), threads=2)
    "Determines the number of steps the MCMC executes"
    nsteps = 5000
    print("Running MCMC...")
    sampler.run_mcmc(p0, nsteps)
    print("Done")

    "Plot the Line-Time"
    plt.clf()
    fig, axes = plt.subplots(4, 1, sharex=True, figsize=(8, 9))
    SaveSamplerChainPlot(axes, 0, sampler.chain[:, :, 0].T, "Flux")
    SaveSamplerChainPlot(axes, 1, sampler.chain[:, :, 1].T+mu_x_correction, "x")
    SaveSamplerChainPlot(axes, 2, sampler.chain[:, :, 2].T+mu_y_correction, "y")
    SaveSamplerChainPlot(axes, 3, sampler.chain[:, :, 3].T, "$\sigma_{x}$")

    "Create Directory for Saving the Plot"
    directory="/Users/edenmolina/PycharmProjects/Quasar/MCMC Plots/%s_Files_%sSteps_%s/" % (SampleName, nsteps, iteration_number)
    if not os.path.exists(directory):
        os.makedirs(directory)
    "Save the Plot"
    fig.savefig(
        "/Users/edenmolina/PycharmProjects/Quasar/MCMC Plots/%s_Files_%sSteps_%s/%s_LineTimePlot_%sSteps_%s.png" % (
        SampleName, nsteps, iteration_number, SampleName, nsteps, iteration_number),
        dpi=264)
    plt.close()

    "Save the Sample Chain Data to Text FIles"
    SaveSampleChainText = False
    if SaveSampleChainText:
        np.savetxt("/Users/edenmolina/PycharmProjects/Quasar/MCMC Plots/%s_Files_%sSteps_%s/%s_SamplerFluxData_%sSteps_%s.txt" % (SampleName, nsteps, iteration_number,SampleName, nsteps, iteration_number), sampler.chain[:, :, 0].T)
        np.savetxt("/Users/edenmolina/PycharmProjects/Quasar/MCMC Plots/%s_Files_%sSteps_%s/%s_SamplerMuXData_%sSteps_%s.txt" % (SampleName, nsteps, iteration_number,SampleName, nsteps, iteration_number), sampler.chain[:, :, 1].T)
        np.savetxt("/Users/edenmolina/PycharmProjects/Quasar/MCMC Plots/%s_Files_%sSteps_%s/%s_SamplerMuYData_%sSteps_%s.txt" % (SampleName, nsteps, iteration_number,SampleName, nsteps, iteration_number), sampler.chain[:, :, 2].T)
        np.savetxt("/Users/edenmolina/PycharmProjects/Quasar/MCMC Plots/%s_Files_%sSteps_%s/%s_SamplerSigXData_%sSteps_%s.txt" % (SampleName, nsteps, iteration_number,SampleName, nsteps, iteration_number), sampler.chain[:, :, 3].T)

    burnin = 500
    samples = sampler.chain[:, burnin:, :].reshape((-1, ndim))
    np.savetxt("/Users/edenmolina/Desktop/SamplerJanuary22.csv", samples)
    "Plotting the corner plot"
    samples.T[1] = samples.T[1] + mu_x_correction
    samples.T[2] = samples.T[2] + mu_y_correction
    fig = corner.corner(samples, labels=["$Flux$", "$X$", "$Y$", "$\sigma_x$"])
    fig.savefig("/Users/edenmolina/PycharmProjects/Quasar/MCMC Plots/%s_Files_%sSteps_%s/%s_Plot_%s_%sSteps.png" % (SampleName, nsteps, iteration_number,SampleName, nsteps, iteration_number), dpi=264)
    plt.close()

    "For Printing Out the Progress"
    width = 10
    for i, result in enumerate(sampler.sample(p0, iterations=nsteps)):
        n = int((width + 1) * float(i) / nsteps)
        sys.stdout.write("\r[{0}{1}]".format('#' * n, ' ' * (width - n)))
    sys.stdout.write("\n")

    "For Printing Out The Medians of the 6 Distribution Functions"
    Names = ['Flux: ', "X: ", "Y: ", "Sig_X: "]
    MedianArray = np.zeros(5)
    for i in range(4):
        MedianArray[i] = np.median(samples.T[i])
        print Names[i], MedianArray[i]
        MedianArray[4] = np.mean(sampler.acceptance_fraction)*100

    #Print out the mean acceptance fraction of the MCMC
    print("Mean acceptance fraction: {0:.3f}"
          .format(np.mean(sampler.acceptance_fraction)*100))

    "Export the Data as CSV"
    df = pd.DataFrame(data = MedianArray, index=['flux', 'mu_x', 'mu_y', 'sigma_x', 'Acceptance Rate'])
    df.to_csv("/Users/edenmolina/PycharmProjects/Quasar/MCMC Plots/%s_Files_%sSteps_%s/%s_Data_%sSteps_%s.csv" % (SampleName, nsteps, iteration_number,SampleName, nsteps, iteration_number))

    "Return the median flux"
    #Use when getting the flux of a 10x10 box
    flux = np.median(samples.T[0])
    return flux
