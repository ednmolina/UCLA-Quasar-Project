# some python imports
__author__ = 'edenm'
import matplotlib.pyplot as plt
import pickle
import cPickle as pickle
from lenstronomy.Plots.output_plots import LensModelPlot
import numpy as np
import pandas as pd
import lenstronomy.Plots.output_plots as out_plot
import scipy.integrate as integrate
import scipy.special as sp

def getPickledData(Path, LensName):
    pick_path = Path

    lens_result = pickle.load(open("%s/lens_result_%s.pickle" % (pick_path, LensName), 'rb'))
    source_result = pickle.load(open("%s/source_result_%s.pickle" % (pick_path, LensName), 'rb'))
    lens_light_result = pickle.load(open("%s/lens_light_result_%s.pickle" % (pick_path, LensName), 'rb'))

    return lens_result, source_result, lens_light_result


def integrand(intensity, r_sersic, n_sersic, r):
    return integrate.quad(2 * np.pi * intensity * r * np.exp( -(r / r_sersic) ^ (1 / n_sersic)), [0, np.inf])

def expintegral(r):
    return quad(integrand, 0, np.inf, args=(r))[0]

def getFlux(Path, LensName):
    "Get the pickled data"
    lens_result, source_result, lens_light_result = getPickledData(
        Path, LensName)

    #define the function of circular sersic
    integrand = lambda r: 2 * np.pi * lens_light_result[0]['amp'] * r * np.exp( -(r / lens_light_result[0]['R_sersic']) ** (1 / lens_light_result[0]['n_sersic']))
    sersic_integral = integrate.quad(integrand, 0, np.inf)

    zero_point = 24.74
    print "Unconverted Flux: ", sersic_integral[0]
    #print "Converted Flux: ", -2.5*np.log10(sersic_integral[0]/zero_point), "\n"
    print -2.5*np.log10(sersic_integral[0]/zero_point)
    #converted_flux = -2.5*(np.log10(sersic_integral/24.74))
    #print converted_flux

    #print 'Lenstronomy Amplitude', lens_light_result[0]['amp']
    #print 'Converted Flux', -2.5 * np.log(lens_light_result[0]['amp']/24.74)
    #print integrate.quad(2 * np.pi* lens_light_result[0]['amp'] * r * np.exp(-(r/lens_light_result[0]['r_sersic'])^(1/lens_light_result[0]['n_sersic'])), [0, np.inf])
    # print lens_result
    # print lens_light_result
    # print ps_result


config_file = pd.read_excel("/Users/edenmolina/PycharmProjects/Quasar/Lenstronomy/LenstronomyConfig.xlsx", sheetname="Master")
Names = config_file['Name']
Paths = config_file['Path']
RunBool = config_file['RunBool']

for i in range(len(Names)):
    if RunBool[i] == 1:
        print "Plotting %s"%Names[i]
        getFlux("/Users/edenmolina/PycharmProjects/Quasar/Lenstronomy", Names[i])
    else:
        print "Skipping %s" % Names[i], "\n"