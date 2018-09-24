# some python imports
__author__ = 'edenm'
import matplotlib.pyplot as plt
import pickle
import cPickle as pickle
from lenstronomy.Plots.output_plots import LensModelPlot
import numpy as np
import pandas as pd
import lenstronomy.Plots.output_plots as out_plot

def getPickledData(Path, LensName):
    pick_path = Path
    lens_result = pickle.load(open("%s/lens_result_%s.pickle" %(pick_path, LensName), 'rb'))
    source_result = pickle.load(open("%s/source_result_%s.pickle" %(pick_path, LensName), 'rb'))
    lens_light_result = pickle.load(open("%s/lens_light_result_%s.pickle" %(pick_path, LensName), 'rb'))
    ps_result = pickle.load(open("%s/ps_result_%s.pickle"%(pick_path, LensName), 'rb'))
    cosmo_result = pickle.load(open("%s/cosmo_result_%s.pickle" %(pick_path, LensName), 'rb'))
    chain_list = pickle.load(open("%s/chain_list_%s.pickle" %(pick_path, LensName), 'rb'))
    param_list = pickle.load(open("%s/param_list_%s.pickle" %(pick_path, LensName), 'rb'))
    samples_mcmc = pickle.load(open("%s/samples_mcmc_%s.pickle" %(pick_path, LensName), 'rb'))
    param_mcmc = pickle.load(open("%s/param_mcmc_%s.pickle" %(pick_path, LensName), 'rb'))
    dist_mcmc = pickle.load(open("%s/dist_mcmc_%s.pickle" %(pick_path, LensName), 'rb'))


    kwargs_data = pickle.load(open("%s/kwargs_data_%s.pickle" %(pick_path, LensName), 'rb'))
    kwargs_psf = pickle.load(open("%s/kwargs_psf_%s.pickle" %(pick_path, LensName), 'rb'))
    kwargs_numerics = pickle.load(open("%s/kwargs_numerics_%s.pickle" %(pick_path, LensName), 'rb'))
    kwargs_model = pickle.load(open("%s/kwargs_model_%s.pickle" %(pick_path, LensName), 'rb'))

    "Import the updates"
    multi_band_list_out = pickle.load(open("%s/multi_band_list_out%s.pickle" %(pick_path, LensName), 'rb'))
    kwargs_params_out = pickle.load(open("%s/kwargs_params_out%s.pickle" %(pick_path, LensName), 'rb'))
    kwargs_psf_out = pickle.load(open("%s/kwargs_psf_out%s.pickle" %(pick_path, LensName), 'rb'))



    return lens_result, source_result, lens_light_result, ps_result, cosmo_result, chain_list, param_list, samples_mcmc, param_mcmc, dist_mcmc, kwargs_data, kwargs_psf, kwargs_numerics, kwargs_model, multi_band_list_out, kwargs_params_out, kwargs_psf_out

def plotPSF_Iteration(kwargs_psf_out, Path, LensName):
    f, axes = out_plot.psf_iteration_compare(kwargs_psf_out, vmin=-6)
    f.savefig("%s/%sPSF_Iteration.pdf" % (Path, LensName))


def getPlots(Path, LensName):
    "Get the pickled data"

    lens_result, source_result, lens_light_result, ps_result, cosmo_result, chain_list, param_list, samples_mcmc, param_mcmc, dist_mcmc, kwargs_data, kwargs_psf, kwargs_numerics, kwargs_model, multi_band_list_out, kwargs_params_out, kwargs_psf_out = getPickledData(Path, LensName)
    "Plotting the Data"
    lensPlot = LensModelPlot(kwargs_data, kwargs_psf, kwargs_numerics, kwargs_model, lens_result, source_result,
                             lens_light_result, ps_result, arrow_size=0.02, cmap_string="gist_heat", high_res=5)

    f, axes = plt.subplots(2, 3, figsize=(16, 8), sharex=False, sharey=False)

    lensPlot.data_plot(ax=axes[0, 0])
    lensPlot.model_plot(ax=axes[0, 1])
    lensPlot.normalized_residual_plot(ax=axes[0, 2], v_min=-6, v_max=6)
    lensPlot.source_plot(ax=axes[1, 0], convolution=False, deltaPix_source=0.01, numPix=100)
    lensPlot.convergence_plot(ax=axes[1, 1], v_max=1)
    lensPlot.magnification_plot(ax=axes[1, 2])
    f.tight_layout()
    f.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0., hspace=0.05)
    plt.savefig("%s/%sPlot1.pdf" %(Path, LensName))
    #plt.show()

    f, axes = plt.subplots(2, 3, figsize=(16, 8), sharex=False, sharey=False)

    lensPlot.decomposition_plot(ax=axes[0, 0], text='Lens light', lens_light_add=True, unconvolved=True)
    lensPlot.decomposition_plot(ax=axes[1, 0], text='Lens light convolved', lens_light_add=True)
    lensPlot.decomposition_plot(ax=axes[0, 1], text='Source light', source_add=True, unconvolved=True)
    lensPlot.decomposition_plot(ax=axes[1, 1], text='Source light convolved', source_add=True)
    lensPlot.decomposition_plot(ax=axes[0, 2], text='All components', source_add=True, lens_light_add=True,
                                unconvolved=True)
    lensPlot.decomposition_plot(ax=axes[1, 2], text='All components convolved', source_add=True, lens_light_add=True,
                                point_source_add=True)
    f.tight_layout()
    f.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0., hspace=0.05)
    plt.savefig("%s/%sPlot2.pdf" %(Path, LensName))
    #plt.show()
    print "This is the lens result", "\n", lens_result #To access a particular value must index the array lens_result[0]
    print "This is the source result", "\n", source_result
    print "This is the lens light result", "\n", lens_light_result
    print "This is the point source result", "\n", ps_result

    import lenstronomy.Plots.output_plots as out_plot
    plt.close()
    for i in range(len(chain_list)):
        if len(param_list[i]) > 0:
            f2, axes2 = out_plot.plot_chain(chain_list[i], param_list[i])
            plt.minorticks_on()
            plt.minorticks_on()
            plt.grid(which="major", alpha=.5, axis='x')
            plt.grid(which="minor", alpha=.5, axis='x', linestyle=":")
            plt.grid(which="major", alpha=.5, axis='y')
            plt.grid(which="minor", alpha=.5, axis='y', linestyle=":")
        f2.savefig("%s/%sPlotParams.pdf" % (Path, LensName))

    import lenstronomy.Util.util as util
    data = lensPlot._data[0:50, 0:50]
    print "Median", np.median(data)
    print "Width", np.std(data)
    data1d = util.image2array(data)
    #plt.hist(data1d, bins=np.linspace(-.2, 0.2, 100))
    #plt.show()
    print "\n"

    "Plot the PSF Iteration"
    plotPSF_Iteration(kwargs_psf_out, Path, LensName)

    "Plot the data and its comonents"
    f, axes = plt.subplots(2, 3, figsize=(16, 8), sharex=False, sharey=False)

    lensPlot.subtract_from_data_plot(ax=axes[0, 0], text='Data')
    lensPlot.subtract_from_data_plot(ax=axes[0, 1], text='Data - Point Source', point_source_add=True)
    lensPlot.subtract_from_data_plot(ax=axes[0, 2], text='Data - Lens Light', lens_light_add=True)
    lensPlot.subtract_from_data_plot(ax=axes[1, 0], text='Data - Source Light', source_add=True)
    lensPlot.subtract_from_data_plot(ax=axes[1, 1], text='Data - Source Light - Point Source', source_add=True,
                                     point_source_add=True)
    lensPlot.subtract_from_data_plot(ax=axes[1, 2], text='Data - Lens Light - Point Source', lens_light_add=True,
                                     point_source_add=True)
    f.tight_layout()
    f.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0., hspace=0.05)
    f.savefig("%s/%sDataComponents.pdf" % (Path, LensName))

config_file = pd.read_excel("/Users/edenmolina/PycharmProjects/Quasar/Lenstronomy/LenstronomyConfig.xlsx", sheetname="Master")
Names = config_file['Name']
Paths = config_file['Path']
RunBool = config_file['RunBool']

for i in range(len(Names)):
    if RunBool[i] == 1:
        print "Plotting %s"%Names[i]
        getPlots("/Users/edenmolina/PycharmProjects/Quasar/Lenstronomy", Names[i])
    else:
        print "Skipping %s" % Names[i], "\n"
