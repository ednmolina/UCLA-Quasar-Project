"""
This is the current workflow that will import a config file that will tell this script how to reconstruct images of lensed quasars
compatible with lenstronomy 0.4.0
"""

# some python imports
__author__ = 'edenm'
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
import time
from mpl_toolkits.axes_grid1 import make_axes_locatable
from lenstronomy.ImSim.image_model import ImageModel
from lenstronomy.Data.imaging_data import Data
from lenstronomy.Data.psf import PSF
from lenstronomy.PointSource.point_source import PointSource
from lenstronomy.LensModel.lens_model import LensModel
from lenstronomy.LensModel.Solver.lens_equation_solver import LensEquationSolver
from lenstronomy.LightModel.light_model import LightModel
# from lenstronomy.Workflow.parameters import Param # Works with version 0.4.0
from lenstronomy.Sampling.parameters import Param # Not currently working with verion 0.4.0
import lenstronomy.Util.util as util
import astropy.io.fits as fits
import seaborn as sns
import lenstronomy.Util.kernel_util as kernel_util
from lenstronomy.Workflow.fitting_sequence import FittingSequence
from lenstronomy.Plots.output_plots import LensModelPlot
import lenstronomy.Util.image_util as im_util
import pickle
import pandas as pd
import ast
import lenstronomy.Plots.output_plots as out_plot


def replaceValuesRADEC(dict, ra_lens, dec_lens):
    for i in dict.values():
        if i == 'ra_lens':
            dict['center_x'] = ra_lens
        if i == 'ra_lens':
            dict['center_y'] = dec_lens
    return dict

def getSquareCutout(source_cood, image, crop_factor):
    x, y = source_cood[0], source_cood[1]
    crop_factor = crop_factor
    cropped_image =  image[int(x-crop_factor):int(x+crop_factor),int(y-crop_factor): int(y+crop_factor)]
    return cropped_image

def rotate_stack_PSF(kernel_point_source):
    kernel_point_source_90 = im_util.rotateImage(kernel_point_source, 90)
    kernel_point_source_180 = im_util.rotateImage(kernel_point_source, 180)
    kernel_point_source_270 = im_util.rotateImage(kernel_point_source, 270)

    mean_kernel_point_source = (kernel_point_source+kernel_point_source_90+kernel_point_source_180+kernel_point_source_270)/(4.)
    return mean_kernel_point_source

def getImage(file_path, exp_time, bkg_rms, bkg_mean, crop_factor, lens_org):
    file_path = file_path

    hdulist = fits.open(file_path)
    image = hdulist[0].data
    #image[image < 0.] = 1e-10
    x_image, y_image = lens_org[0], lens_org[1]

    cropped_image = getSquareCutout([x_image, y_image], image, crop_factor)

    plt.imshow(np.log10(cropped_image),  # use log10 as the scale than let's you compare change in magnitude
               origin='lower',
               vmin=-.9,
               vmax=2,
               cmap=sns.cubehelix_palette(start=0.5, rot=-1.5, gamma=1., hue=1.,
                                          # l/ight=.8,
                                          dark=0.,
                                          reverse=True, as_cmap=True))
    plt.colorbar()
    plt.title("Cropped Image")
    plt.grid(False)
    #plt.show()

    numPix = len(cropped_image)
    print ("THIS IS NUMPIX", numPix)
    deltaPix = 0.009942
    x_grid, y_grid, ra_at_xy_0, dec_at_xy_0, x_at_radec_0, y_at_radec_0, Mpix2coord, Mcoord2pix = util.make_grid_with_coordtransform(
        numPix=numPix, deltapix=deltaPix, subgrid_res=1)
    # mask (1= model this pixel, 0= leave blanck)
    # replace with wht_frame (in fits file data[1]) and CCD_gain (header)
    exposure_time = exp_time

    exposure_map = np.ones((numPix, numPix)) * exposure_time  # individual exposure time/weight per pixel
    background_rms = bkg_rms  # estimate from empty patch of sky
    mean_bkg = bkg_mean # mean background flux; X^2 = .53

    kwargs_data = {
        'background_rms': background_rms,
        'exposure_map': exposure_map
        , 'ra_at_xy_0': ra_at_xy_0, 'dec_at_xy_0': dec_at_xy_0, 'transform_pix2angle': Mpix2coord
        , 'image_data': cropped_image - mean_bkg
    }
    data_class = Data(kwargs_data)
    return data_class, deltaPix, kwargs_data, cropped_image, image

def getPSF(data_class, deltaPix, kwargs_data, x_quasar, y_quasar, lens_cropped, kernel_size, Path, LensName, PSF_Choice, PSF_Int_Kernel, PSF_Rot):
    x_quasar = x_quasar
    y_quasar = y_quasar

    ra_quasar, dec_quasar = data_class.map_pix2coord(x_quasar, y_quasar)
    theta_e_est = ((((x_quasar[1] - x_quasar[0]) ** 2 + (y_quasar[1] - y_quasar[0]) ** 2) ** .5) / 2) * deltaPix

    x_lens = lens_cropped[0]
    y_lens = lens_cropped[1]

    ra_lens, dec_lens = data_class.map_pix2coord(x_lens, y_lens)
    print ("This is the converted RA and Dec of Lens Line 112: ", ra_lens, dec_lens)

    print "Theta_E Estimate: %s" % theta_e_est

    kernel_size = kernel_size
    print "This is the kernel size for PSF before padding with 0's", kernel_size

    kernel_point_source = kernel_util.cutout_source(x_quasar[PSF_Choice], y_quasar[PSF_Choice], kwargs_data['image_data'],
                                                    kernelsize=kernel_size, shift=True)
    kernel_point_source /= np.sum(kernel_point_source)

    if PSF_Rot == 1:
        kernel_point_source = rotate_stack_PSF(kernel_point_source)

    num_large = int(kernel_size+kernel_size*.5)
    if num_large%2 == 0:
        num_large = num_large + 1
    else:
        num_large = num_large

    kernel_point_source_new = np.zeros((num_large, num_large))


    num_small = len(kernel_point_source)
    print "THIS IS NUM SMALL", num_small
    x_min = (num_large - num_small) / 2
    x_max = x_min + num_small

    kernel_point_source_new[x_min:x_max, x_min:x_max] = kernel_point_source

    psf_error_map = np.ones_like(kernel_point_source_new) * PSF_Int_Kernel #Edita

    np.savetxt("/Users/edenmolina/Desktop/PSF_ErrorMap.txt", psf_error_map)

    kwargs_psf = {'psf_type': 'PIXEL', 'kernel_point_source': kernel_point_source_new,
                  'kernel_point_source_init': kernel_point_source_new,
                  'psf_error_map': psf_error_map}

    "For plotting the PSF"
    plt.close()
    plt.imshow(np.log10(kernel_point_source), origin='lower')
    plt.title("PSF Error Map")
    plt.savefig("%s/%s_Plot0_PSF.pdf" % (Path, LensName))

    plt.close()

    return kwargs_psf, ra_quasar, dec_quasar, ra_lens, dec_lens, theta_e_est

def getPSF_MIN(data_class, deltaPix, kwargs_data, x_quasar, y_quasar, lens_cropped, kernel_size, Path, LensName, PSF_Choice, PSF_Int_Kernel, PSF_Rot):
    x_quasar = x_quasar
    y_quasar = y_quasar

    ra_quasar, dec_quasar = data_class.map_pix2coord(x_quasar, y_quasar)
    theta_e_est = ((((x_quasar[1] - x_quasar[0]) ** 2 + (y_quasar[1] - y_quasar[0]) ** 2) ** .5) / 2) * deltaPix

    x_lens = lens_cropped[0]
    y_lens = lens_cropped[1]

    ra_lens, dec_lens = data_class.map_pix2coord(x_lens, y_lens)
    print ("This is the converted RA and Dec of Lens Line 175: ", ra_lens, dec_lens)

    print "Theta_E Estimate: %s" % theta_e_est

    kernel_size = kernel_size
    print "This is the kernel size for PSF before padding with 0's", kernel_size

    kernel_point_source = kernel_util.cutout_source(x_quasar[PSF_Choice], y_quasar[PSF_Choice], kwargs_data['image_data'],
                                                    kernelsize=kernel_size, shift=True)

    #Get the PSF choice here; essentially chooses which of the two quasars to use as a PSF
    if PSF_Choice == 1:
        PSF_Choice = 0
        kernel_point_source_other = kernel_util.cutout_source(x_quasar[PSF_Choice], y_quasar[PSF_Choice],
                                                        kwargs_data['image_data'],
                                                        kernelsize=kernel_size, shift=True)
    if PSF_Choice == 0:
        PSF_Choice = 1
        kernel_point_source_other = kernel_util.cutout_source(x_quasar[PSF_Choice], y_quasar[PSF_Choice],
                                                              kwargs_data['image_data'],
                                                              kernelsize=kernel_size, shift=True)

    kernel_point_source /= np.sum(kernel_point_source)
    kernel_point_source_other /= np.sum(kernel_point_source_other)

    if PSF_Rot == 1:
        kernel_point_source = rotate_stack_PSF(kernel_point_source)
        kernel_point_source_other = rotate_stack_PSF(kernel_point_source_other)

    num_large = int(kernel_size+kernel_size*.5)
    if num_large%2 == 0:
        num_large = num_large + 1
    else:
        num_large = num_large

    kernel_point_source_new = np.zeros((num_large, num_large))
    kernel_point_source_other_new = np.zeros((num_large, num_large))

    num_small = len(kernel_point_source)
    print "THIS IS NUM SMALL", num_small
    x_min = (num_large - num_small) / 2
    x_max = x_min + num_small

    #Crop the images
    kernel_point_source_new[x_min:x_max, x_min:x_max] = kernel_point_source
    kernel_point_source_other_new[x_min:x_max, x_min:x_max] = kernel_point_source_other

    #Normalize the two PSF cutouts
    kernel_point_source_new = kernel_point_source_new/np.sum(kernel_point_source_new)
    kernel_point_source_other_new = kernel_point_source_other_new / np.sum(kernel_point_source_other_new)

    "NEED TO FIND THE ELEMENTAL MIN OF TWO PSF's"
    kernel_point_source_final = np.minimum(kernel_point_source_new, kernel_point_source_other_new)
    psf_error_map = np.ones_like(kernel_point_source_final) * PSF_Int_Kernel #Edita

    # np.savetxt("/Users/edenmolina/Desktop/PSF_ErrorMap.txt", psf_error_map)

    kwargs_psf = {'psf_type': 'PIXEL', 'kernel_point_source': kernel_point_source_final,
                  'kernel_point_source_init': kernel_point_source_final
                   ,'psf_error_map': psf_error_map
                  }

    "For plotting the PSF"
    plt.close()
    plt.imshow(np.log10(kernel_point_source_final), origin='lower')
    plt.title("PSF Error Map")
    plt.savefig("%s/%s_Plot0_PSF.pdf" % (Path, LensName))
    #plt.show()
    plt.close()

    return kwargs_psf, ra_quasar, dec_quasar, ra_lens, dec_lens, theta_e_est

def setParameters(theta_e_est, ra_lens, dec_lens, ra_quasar, dec_quasar, lens_params_dict, source_params_dict, lens_light_params_dict, ps_params_dict):

    "Edit the dictionary to replace the placeholders for the strings: ra/dec_lens"
    lens_params_dict = replaceValuesRADEC(lens_params_dict, ra_lens, dec_lens)
    source_params_dict = replaceValuesRADEC(source_params_dict, ra_lens, dec_lens)

    # initial guess of non-linear parameters, we chose different starting parameters than the truth #
    kwargs_lens_init = [
        {'theta_E': theta_e_est, 'e1': 0., 'e2': 0., 'gamma': 2., 'center_x': ra_lens, 'center_y': dec_lens},
        {'e1': 0., 'e2': 0.}]
    kwargs_source_init = [
        {'R_sersic': 0.02, 'n_sersic': 2., 'e1': 0.2, 'e2': 0.1, 'center_x': ra_lens, 'center_y': dec_lens}]
    kwargs_lens_light_init = [
        {'R_sersic': 0.5, 'n_sersic': 2., 'e1': 0.1, 'e2': 0.1, 'center_x': ra_lens, 'center_y': dec_lens}]
    kwargs_ps_init = [{'ra_image': ra_quasar+0.01, 'dec_image': dec_quasar-0.01}]

    # initial spread in parameter estimation #
    kwargs_lens_sigma = [{'theta_E': 0.1, 'e1': 0.1, 'e2': 0.1, 'gamma': .1, 'center_x': 0.1, 'center_y': 0.1},
                         {'e1': 0.1, 'e2': 0.1}]
    kwargs_source_sigma = [{'R_sersic': 0.2, 'n_sersic': .5, 'center_x': .1, 'center_y': 0.2, 'e1': 0.2, 'e2': 0.2}]
    kwargs_lens_light_sigma = [
        {'R_sersic': 0.3, 'n_sersic': 0.5, 'e1': 0.2, 'e2': 0.2, 'center_x': 1, 'center_y': 1}]
    kwargs_ps_sigma = [{'ra_image': [0.02] * 2, 'dec_image': [0.02] * 2}]

    # hard bound lower limit in parameter space #
    kwargs_lower_lens = [{'theta_E': .01, 'e1': -0.2, 'e2': -0.2, 'gamma': 1.5, 'center_x': -5., 'center_y': -5},
                         {'e1': -.2, 'e2': -0.2}]
    kwargs_lower_source = [
        {'R_sersic': 0.01, 'n_sersic': 0.5, 'e1': -0.2, 'e2': -0.2, 'center_x': -2., 'center_y': -2.}]
    kwargs_lower_lens_light = [
        {'R_sersic': 0.1, 'n_sersic': 0.5, 'e1': -0.2, 'e2': -0.2, 'center_x': -1., 'center_y': -1.}]
    kwargs_lower_ps = [{'ra_image': -10 * np.ones_like(ra_quasar), 'dec_image': -10. * np.ones_like(dec_quasar)}]

    # hard bound upper limit in parameter space #
    kwargs_upper_lens = [{'theta_E': 2, 'e1': 0.2, 'e2': 0.2, 'gamma': 2.5, 'center_x': 5., 'center_y': 5.},
                         {'e1': 0.2, 'e2': 0.2}]
    kwargs_upper_source = [{'R_sersic': 2., 'n_sersic': 8., 'e1': 0.2, 'e2': 0.2, 'center_x': 2., 'center_y': 2.}]
    kwargs_upper_lens_light = [{'R_sersic': 4, 'n_sersic': 8., 'e1': 0.2, 'e2': 0.2, 'center_x': 2., 'center_y': 2.}]
    kwargs_upper_ps = [{'ra_image': 10 * np.ones_like(ra_quasar), 'dec_image': 10 * np.ones_like(dec_quasar)}]

    "Can fix parameters here"
    lens_params = [kwargs_lens_init, kwargs_lens_sigma, [lens_params_dict, {}], kwargs_lower_lens,
                   kwargs_upper_lens]
    source_params = [kwargs_source_init, kwargs_source_sigma, [source_params_dict ,{}], kwargs_lower_source, kwargs_upper_source]
    lens_light_params = [kwargs_lens_light_init, kwargs_lens_light_sigma, [lens_light_params_dict, {}], kwargs_lower_lens_light,
                         kwargs_upper_lens_light]
    ps_params = [kwargs_ps_init, kwargs_ps_sigma, [ps_params_dict, {}], kwargs_lower_ps, kwargs_upper_ps]

    print ("These are the parameters")
    print (lens_params, "\n", source_params, "\n", lens_light_params,"\n", ps_params)
    print ("These are the parameter inputs (dictionaries). These are also the values which are fixed in the model.")
    print (lens_params_dict, source_params_dict, lens_light_params_dict, ps_params_dict)
    # kwargs_init = [kwargs_lens, kwargs_source, kwargs_lens_light, kwargs_ps]
    kwargs_params = {'lens_model': lens_params,
                     'source_model': source_params,
                     'lens_light_model': lens_light_params,
                     'point_source_model': ps_params}

    return kwargs_params, kwargs_lens_sigma, kwargs_lower_lens, kwargs_upper_lens, kwargs_source_sigma, kwargs_lower_source, kwargs_upper_source, kwargs_lens_light_sigma, kwargs_lower_lens_light, kwargs_upper_lens_light, kwargs_ps_sigma, kwargs_lower_ps, kwargs_upper_ps

"Initialize Lenstronomy"
def lenstronomy_master(LensName, file_path, exp_time, bkg_rms, bkg_mean, x_quasar, y_quasar, crop_factor, lens_org, lens_cropped, kernel_size, lens_params_dict, source_params_dict, lens_light_params_dict, ps_params_dict, n_iterations, n_particles, PSF_Choice, PSF_Int_Kernel, SourceVisible, PSF_Rot, SigmaScale, *mask_path):

    "Importing the image file"
    data_class, deltaPix, kwargs_data, cropped_image, image = getImage(file_path, exp_time, bkg_rms, bkg_mean, crop_factor, lens_org)
    Path = "/Users/edenmolina/PycharmProjects/Quasar/Lenstronomy"

    "Set the Save Destination"
    pick_path = "/Users/edenmolina/PycharmProjects/Quasar/Lenstronomy"

    "Getting the PSF of a lensed Quasar"
    kwargs_psf, ra_quasar, dec_quasar, ra_lens, dec_lens, theta_e_est = getPSF_MIN(data_class, deltaPix, kwargs_data, x_quasar, y_quasar, lens_cropped, kernel_size, Path, LensName, PSF_Choice, PSF_Int_Kernel, PSF_Rot)

    "Lens Model Parameters"
    lens_model_list = ['SIE', 'SHEAR']
    kwargs_shear = {'e1': 0.01, 'e2': 0.01}  # gamma_ext: shear strength, psi_ext: shear angel (in radian)
    kwargs_sie = {'theta_E': theta_e_est, 'center_x': 0, 'center_y': 0}
    kwargs_lens = [kwargs_sie, kwargs_shear]
    lens_model_class = LensModel(lens_model_list=lens_model_list)

    #source_model_list = ['SERSIC']
    "Model a soource if it is seen"
    if SourceVisible == 0:
        source_model_list = []
    else:
        source_model_list = ['SERSIC']

    lens_light_model_list = ['SERSIC_ELLIPSE']  # , SERSIC_ELLIPSE]
    point_source_list = ['LENSED_POSITION']

    kwargs_model = {'lens_model_list': lens_model_list,
                   'source_light_model_list': source_model_list,
                   'lens_light_model_list': lens_light_model_list,
                   'point_source_model_list': point_source_list,
                   'fixed_magnification_list': [False],
                   }

    "Checks to see if there is a mask for masking out parts of the images such as extra stars"
    if mask_path==(0,):
        print "No Mask"
        #kwargs_numerics = {'subgrid_res': 1, 'psf_subgrid': False, 'mask': np.ones_like(cropped_image)}
        kwargs_numerics = {'subgrid_res': 1, 'psf_subgrid': False,
                           'mask': np.ones_like(cropped_image),
                           'psf_keep_error_map': True,
                           'point_source_subgrid': 1}
        plt.imshow(np.ones_like(cropped_image), origin = 'lower')
        #plt.show()
        multi_band_list = [[kwargs_data, kwargs_psf, kwargs_numerics]]
    else:
        print "Mask found"
        hdulist = fits.open(str(mask_path[0]))
        mask = hdulist[0].data

        cropped_mask = getSquareCutout(lens_org, mask, int(crop_factor))
        cropped_mask = 1-cropped_mask
        print cropped_mask
        # plt.close()
        plt.imshow(cropped_mask, origin='lower')
        #plt.show()
        kwargs_numerics = {'subgrid_res': 1, 'psf_subgrid': False,
                           'mask': cropped_mask,
                           'psf_keep_error_map': True,
                           'point_source_subgrid': 1}
        multi_band_list = [[kwargs_data, kwargs_psf, kwargs_numerics]]

    num_source_model = len(source_model_list)

    if SourceVisible == 0:
        kwargs_constraints = {
           'num_point_source_list': [2],
           'additional_images_list': [False],
           'joint_lens_with_light': [[0, 0, ['center_x', 'center_y']]],
           'joint_source_with_point_source': [],  # [[0, 0]],
           # 'fix_to_point_source_list': [True] * num_source_model,
           'image_plane_source_list': [False] * num_source_model,
           #'joint_lens_light_with_lens_light': [[0, 1, ['center_x', 'center_y']]]
           'solver': True,
           'solver_type': 'THETA_E_PHI',  # 'PROFILE', 'PROFILE_SHEAR', 'ELLIPSE', 'CENTER'

           }

        pickle.dump(kwargs_constraints, open("%s/kwargs_constraintsSI.pickle" % (pick_path), 'wb'))

    else:

        kwargs_constraints = {
           'num_point_source_list': [2],
           'additional_images_list': [False],
           'joint_lens_with_light': [[0, 0, ['center_x', 'center_y']]],
           'joint_source_with_point_source': [[0, 0]],
           #'fix_to_point_source_list': [True] * num_source_model,
           'image_plane_source_list': [False] * num_source_model,
           'joint_lens_light_with_lens_light': [[0, 0, ['center_x', 'center_y']]],
           'solver': True,
           'solver_type': 'THETA_E_PHI',  # 'PROFILE', 'PROFILE_SHEAR', 'ELLIPSE', 'CENTER'
           }

        pickle.dump(kwargs_constraints, open("%s/kwargs_constraintsSV.pickle" % (pick_path), 'wb'))

    kwargs_likelihood = {'check_bounds': True,
                        'force_no_add_image': False,
                        'source_marg': False,
                        'point_source_likelihood': False,
                        'position_uncertainty': 0.004,
                        'check_solver': True,
                        'solver_tolerance': 0.001,
                        #'force_no_add_image': True #Ensure that quads don't occur
                        'restrict_image_number': True, #Also to ensure that quads are not modeled
                        'check_positive_flux': True # Check if getting negative flux
                        }

    "Initialization of parameters for parameter space"

    kwargs_params, kwargs_lens_sigma, kwargs_lower_lens, kwargs_upper_lens, kwargs_source_sigma, kwargs_lower_source, kwargs_upper_source, kwargs_lens_light_sigma, kwargs_lower_lens_light, kwargs_upper_lens_light, kwargs_ps_sigma, kwargs_lower_ps, kwargs_upper_ps = setParameters(theta_e_est, ra_lens, dec_lens, ra_quasar, dec_quasar, lens_params_dict, source_params_dict, lens_light_params_dict, ps_params_dict)

    """
            ******Initialize the PSO and MCMC Parameters******
    """
    "First fitting sequence"
    fitting_seq = FittingSequence(multi_band_list, kwargs_model, kwargs_constraints, kwargs_likelihood, kwargs_params)

    # PSF Reconstruction Occurs here
    psf_iterations = 75
    fitting_kwargs_list = [
       {'fitting_routine': 'PSO', 'mpi': False, 'sigma_scale': 1., 'n_particles': n_particles,
        'n_iterations': n_iterations},
       {'fitting_routine': 'psf_iteration', 'psf_iter_num': psf_iterations, 'psf_iter_factor': 0.1},
   ]

    lens_result, source_result, lens_light_result, ps_result, cosmo_result, chain_list, param_list, samples_mcmc, param_mcmc, dist_mcmc = fitting_seq.fit_sequence(
       fitting_kwargs_list)

    multi_band_list_out = fitting_seq.multi_band_list
    kwargs_params_out = fitting_seq.kwargs_params
    kwargs_data, kwargs_psf_out, kwargs_numerics = multi_band_list_out[0]

    print "Saving the Chainlist"
    pickle.dump(chain_list, open("%s/chain_list_%s.pickle" % (pick_path, LensName), 'wb'))
    pickle.dump(param_list, open("%s/param_list_%s.pickle" % (pick_path, LensName), 'wb'))

    "Second Fitting Sequence"
    # The second PSF Reconstruction Occurs here to improve from the first
    lens_params = [lens_result, kwargs_lens_sigma, [lens_params_dict, {}], kwargs_lower_lens, kwargs_upper_lens]
    source_params = [source_result, kwargs_source_sigma, [{}], kwargs_lower_source, kwargs_upper_source]
    lens_light_params = [lens_light_result, kwargs_lens_light_sigma, [{}], kwargs_lower_lens_light,
                        kwargs_upper_lens_light]
    ps_params = [ps_result, kwargs_ps_sigma, [{}], kwargs_lower_ps, kwargs_upper_ps]

    kwargs_params_out = {'lens_model': lens_params,
                        'source_model': source_params,
                        'lens_light_model': lens_light_params,
                        'point_source_model': ps_params}


    fitting_seq = FittingSequence(multi_band_list_out, kwargs_model, kwargs_constraints, kwargs_likelihood,
                                 kwargs_params_out)


    """--------Running the MCMC --------"""
    fitting_kwargs_list = [

       {'fitting_routine': 'psf_iteration', 'psf_iter_num': 75, 'psf_iter_factor': 0.1},

       ## Comment from Anowar: Uncomment the following line to do MCMC, first do an MCMC for n_burn=0 and n_run about 500, to see where
       ## chains are burnt in, then you can choose n_burn and n_run wisely. The MCMC chain will be returned in the
       ##  samples_mcmc variable in the line below. Save samples_mcmc in a file as other variables.
       {'fitting_routine': 'MCMC', 'n_burn': 0, 'n_run': 700, 'walkerRatio': 10, 'mpi': False,'sigma_scale': SigmaScale}
   ]
    t0 = time.time()
    lens_result, source_result, lens_light_result, ps_result, cosmo_result, chain_list, param_list, samples_mcmc, param_mcmc, dist_mcmc = fitting_seq.fit_sequence(
       fitting_kwargs_list)
    t1 = time.time()
    print ("It took %s seconds to run the MCMC." %(t1-t0))
    multi_band_list_out = fitting_seq.multi_band_list
    kwargs_params_out = fitting_seq.kwargs_params
    kwargs_data, kwargs_psf_out, kwargs_numerics = multi_band_list_out[0]


    "Saving the outputs"
    print "SOURCE RESULT"
    print source_result
    pickle.dump(kwargs_data, open("%s/kwargs_data_%s.pickle" % (pick_path, LensName), 'wb'))
    pickle.dump(kwargs_psf_out, open("%s/kwargs_psf_out_%s.pickle" % (pick_path, LensName), 'wb'))
    pickle.dump(kwargs_numerics, open("%s/kwargs_numerics_%s.pickle" % (pick_path, LensName), 'wb'))
    pickle.dump(kwargs_model, open("%s/kwargs_model_%s.pickle" % (pick_path, LensName), 'wb'))
    pickle.dump(lens_result, open("%s/lens_result_%s.pickle" % (pick_path, LensName), 'wb'))
    pickle.dump(source_result, open("%s/source_result_%s.pickle" % (pick_path, LensName), 'wb'))
    pickle.dump(lens_light_result, open("%s/lens_light_result_%s.pickle" % (pick_path, LensName), 'wb'))
    pickle.dump(ps_result, open("%s/ps_result_%s.pickle" % (pick_path, LensName), 'wb'))

    # Save the MCMC File
    pickle.dump(samples_mcmc, open("%s/samples_mcmc_result_%s.pickle" % (pick_path, LensName), 'wb'))
    pickle.dump(param_mcmc, open("%s/param_mcmc_result_%s.pickle" % (pick_path, LensName), 'wb'))

    #Save the kwargs params out
    pickle.dump(kwargs_params_out, open("%s/kwargs_params_out_result_%s.pickle" % (pick_path, LensName), 'wb'))


    print "Finished!", "\n"


"IMPORT THE CONFIGURATION FILE AND INITIALIZE THE PARAMETERS FOR THE MODELING"
#lenstronomy_master('0047', '/Users/edenmolina/Documents/Astro/Coadds/2013/0047all.fits', 2340, 0.224522621129, .01650891684, [149, 171], [90, 231], 150, [609, 628], [147, 152], 141, {'gamma': 2, 'e1': 0}, {'n_sersic': 2}, {'n_sersic':2}, {}, 20 , 100)
config_file = pd.read_excel("/Users/edenmolina/PycharmProjects/Quasar/Lenstronomy/LenstronomyConfig.xlsx", sheetname="Master")

Names = config_file['Name']
Paths = config_file['Path']
ExpTime = config_file['ExpTime']
BKG_RMS = config_file['RMS']
BKG_MEAN = config_file['Mean']
Quasar1_x = config_file['Quasar1_x']
Quasar2_x = config_file['Quasar2_x']
Quasar1_y = config_file['Quasar1_y']
Quasar2_y = config_file['Quasar2_y']
Crop_Factor = config_file['Crop_Factor']
LensXOrg = config_file['LensXOrg']
LensYOrg = config_file['LensYOrg']
LensXCrop = config_file['LensXCrop']
LensYCrop = config_file['LensYCrop']

kernel_size = config_file['kernel_size']
lens_params_dict = config_file['Lens_params_dict'].values
source_params_dict = config_file['source_params_dict'].values
lens_light_params_dict = config_file['lens_light_params_dict'].values
ps_params_dict = config_file['ps_params_dict'].values
n_iterations = config_file['n_iterations']
n_particles = config_file['n_particles']
PSF_Choice = config_file['PSF_Choice']
PSF_Int_Kernel = config_file['PSF_Int_kernel']
Mask = config_file['Mask']
SourceVisible = config_file['SourceVisble']
PSF_Rot = config_file['PSF_Rot']
RunBool = config_file['RunBool']
SigmaScale = config_file['SigmaScale']

for i in range(len(Names)):
    if RunBool[i] == 1:
        print "\n", "Working on %s"%Names[i], "\n"

        lenstronomy_master(Names[i], Paths[i], ExpTime[i], BKG_RMS[i], BKG_MEAN[i], [Quasar1_x[i], Quasar2_x[i]],
                           [Quasar1_y[i], Quasar2_y[i]], Crop_Factor[i], [LensXOrg[i], LensYOrg[i]],
                           [LensXCrop[i], LensYCrop[i]], kernel_size[i],ast.literal_eval(lens_params_dict[i]), ast.literal_eval(source_params_dict[i]),
                           ast.literal_eval(lens_light_params_dict[i]), ast.literal_eval(ps_params_dict[i]), n_iterations[i], n_particles[i], PSF_Choice[i], PSF_Int_Kernel[i], SourceVisible[i], PSF_Rot[i], SigmaScale[i], Mask[i])
    else:
        print "Skipping %s" %Names[i], "\n"
