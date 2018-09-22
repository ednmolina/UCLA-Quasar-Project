
"Stable working version"
"""
What's new:
            Now Crops the images and processes a smaller cutout of the original image with the two quasar images nicely cente
            Now also interpoaltes the images to increase precision in the PSF subtraction
This scipt takes in a config file that will perform PSF subtraction on a doubly lensed quasar system
    Takes two quasars, normalizes their flux (determined by the MCMC scrips) and essentially removes them from the image in hopes of revealing a faint galaxy
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as scinterp
import astropy.io.fits as fits
import pylab
import pandas as pd
import sys
import os
import matplotlib as mpl
mpl.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt
import seaborn as sns
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = 'Times New Roman'
plt.rcParams['mathtext.rm'] = 'serif'
plt.rcParams['mathtext.it'] = 'serif:italic'
plt.rcParams['mathtext.bf'] = 'serif:bold'
plt.rcParams['mathtext.fontset'] = 'custom'

sns.set(style='ticks', context='paper', font='Times New Roman', font_scale=2.)
sns.set_style({"xtick.direction": "in", "ytick.direction": "in", "axes.linewidth": 1.0, })

"To Make Cutout of the Quasar/Galaxy for PSF"
def circle_mask(image, object_coor, radius):
    image_data = image.copy()
    nx, ny = image_data.shape
    x, y = np.mgrid[0:nx, 0:ny]
    r = np.sqrt((x-object_coor[1])**2 + (y-object_coor[0])**2)

    for i in range(r.shape[0]):
        for j in range(r.shape[1]):
            if r[i][j] > radius:
                image_data[i][j] = 0
            else:
                image_data[i][j] = image_data[i][j]
    return image_data

def square_mask(image, object_coor, radius):
    image_data = image.copy()
    nx, ny = image_data.shape
    square_x_min = object_coor[0] - .5 * radius
    square_y_min = object_coor[1] - .5 * radius
    square_x_max = object_coor[0] + .5 * radius
    square_y_max = object_coor[1] + .5 * radius

    for x in range(nx):
        for y in range(ny):
            if  square_x_min < x < square_x_max and square_y_min < y < square_y_max:
                image_data[x][y] = image_data[x][y]
            else:
                image_data[x][y] = 0

    return image_data

"Get the proper shift for the two objects based on their relative position"
def getshift(obj1_coord, obj2_coord):
    obj1 = obj1_coord.copy()
    obj2 = obj2_coord.copy()

    x1, y1 = obj1[1], obj1[0]
    x2, y2 = obj2[1], obj2[0]

    x_diff = np.abs(x1 - x2)
    y_diff = np.abs(y1 - y2)
    if x2 > x1:
        x1_shift =  x_diff
        x2_shift = - x_diff

        xobj_min = x1
        xobj_max = x2

        if y1 > y2:
            y1_shift = - y_diff
            y2_shift = y_diff

            yobj_min = y2
            yobj_max = y1

            return x1_shift, x2_shift, y1_shift, y2_shift, int(xobj_min), int(xobj_max), int(yobj_min), int(yobj_max)
        elif y2 > y1:
            y1_shift =  y_diff
            y2_shift = - y_diff

            yobj_min = y1
            yobj_max = y2

            return x1_shift, x2_shift, y1_shift, y2_shift, int(xobj_min), int(xobj_max), int(yobj_min), int(yobj_max)
        elif y2 == y1:
            y1_shift = 0
            y2_shift = 0

            yobj_min = y2
            yobj_max = y1

            return x1_shift, x2_shift, y1_shift, y2_shift, int(xobj_min), int(xobj_max), int(yobj_min), int(yobj_max)

    elif x2 < x1:
        x1_shift = - x_diff
        x2_shift = + x_diff

        xobj_min = x2
        xobj_max = x2
        if y1 > y2:
            y1_shift = - y_diff
            y2_shift = y_diff

            yobj_min = y2
            yobj_max = y1

            return x1_shift, x2_shift, y1_shift, y2_shift, int(xobj_min), int(xobj_max), int(yobj_min), int(yobj_max)
        elif y2 > y1:
            y1_shift =  y_diff
            y2_shift = - y_diff

            yobj_min = y1
            yobj_max = y2

            return x1_shift, x2_shift, y1_shift, y2_shift, int(xobj_min), int(xobj_max), int(yobj_min), int(yobj_max)
        elif y2 == y1:
            y1_shift = 0
            y2_shift = 0

            yobj_min = y2
            yobj_max = y1

            return x1_shift, x2_shift, y1_shift, y2_shift, int(xobj_min), int(xobj_max), int(yobj_min), int(yobj_max)
    elif x2 == x1:
        x1_shift = 0
        x2_shift = 0

        xobj_min = x2
        xobj_max = x2

        if y1 > y2:
            y1_shift = - y_diff
            y2_shift = y_diff

            yobj_min = y2
            yobj_max = y1

            return x1_shift, x2_shift, y1_shift, y2_shift, int(xobj_min), int(xobj_max), int(yobj_min), int(yobj_max)
        elif y2 > y1:
            y1_shift =  y_diff
            y2_shift = - y_diff

            yobj_min = y1
            yobj_max = y2

            return x1_shift, x2_shift, y1_shift, y2_shift, int(xobj_min), int(xobj_max), int(yobj_min), int(yobj_max)
        elif y2 == y1:
            y1_shift = 0
            y2_shift = 0

            yobj_min = y2
            yobj_max = y1

            return x1_shift, x2_shift, y1_shift, y2_shift, int(xobj_min), int(xobj_max), int(yobj_min), int(yobj_max)

"Get the amount to crop the image by"
def getZoom(obj1_coord, obj2_coord, image):
    #Set the coordinates
    obj1_coordcopy = np.copy(obj1_coord)
    obj2_coordcopy = np.copy(obj2_coord)

    x1, y1 = obj1_coordcopy[1], obj1_coordcopy[0]
    x2, y2 = obj2_coordcopy[1], obj2_coordcopy[0]

    im_x, im_y = image.shape[0], image.shape[1]

    minlist = list()
    "Perform comparisons to determine which of the 4 distances from the two images is the smallest"
    if x1 < x2:
        d_x1 = x1
        d_x2 = x2 - im_x
        minlist.append(min(d_x1, d_x2))
    else:
        d_x1 = x1- im_x
        d_x2 = x2
        minlist.append(min(d_x1, d_x2))
    if y1 < y2:
        d_y1 = x1
        d_y2 = x2 - im_y
        minlist.append(min(d_y1, d_y2))
    else:
        d_y1 = y1- im_y
        d_y2 = y2
        minlist.append(min(d_y1, d_y2))
    return np.abs(min(minlist))

def getImageSeparation(obj1_coord, obj2_coord):
    distance = np.sqrt((obj1_coord[1]-obj2_coord[1])**2 + (obj1_coord[0]-obj2_coord[0])**2)
    return distance

def getMinimizedQuadratic(x_Norm_space, Chi_Sq_Array):
    "Get the constnats for the quadratic fit"
    #constants[0] is quadratic term, constants[1] is linear term, constnats[2] is constant term
    constants = np.polyfit(x_Norm_space, Chi_Sq_Array, 2)

    x_min = -constants[1]/(2*constants[0])
    y_min = constants[0] * (x_min**2) + constants[1] * x_min + constants[2]
    return x_min, y_min, constants

def plotChiSq(x_Norm_space, constants, x_min, y_min, objectname, Chi_Sq_Array):
    "Plot Chi Squared as a function of the normalization constant"
    x_plot = np.linspace(x_Norm_space[0], x_Norm_space[-1], 100)
    Quadratic_Eqn = constants[0]*x_plot**2 + constants[1]*x_plot + constants[2]
    plt.close()
    plt.plot(x_plot, Quadratic_Eqn)  # Plots the various narmalization constants possible
    plt.scatter(x_min, y_min, label="Minimum at %s" % round(x_min, 2))
    plt.scatter(x_Norm_space, Chi_Sq_Array)
    plt.xlabel("Normalization Constant")
    plt.ylabel("$\chi^2$")
    plt.legend()
    plt.savefig("/Users/edenmolina/Desktop/New PSF OutPut/%s_ChiSq.pdf" % objectname)
    plt.close()

def getChi_Sq(data_cutout, im_cutout, im1overim2, Object_coord, xy_min, objectname):
    xmin = xy_min[0]
    ymin = xy_min[1]

    # The new coordiantes of the objects with respect to the higher resolution image
    x_im = (Object_coord[0] - (xmin))
    y_im = (Object_coord[1] - (ymin))

    # Coordinates from which we will perform the Chi Sq Minimization
    zoom = 5
    x_min = int(x_im - zoom)
    x_max = int(x_im + zoom)
    y_min = int(y_im - zoom)
    y_max = int(y_im + zoom)

    cutout_data = data_cutout[y_min:y_max, x_min:x_max]
    cutout_im_cutout = im_cutout[y_min:y_max, x_min:x_max]

    "THIS IS TEMP"
    plt.close()
    plt.title("Image cutout for %s"%objectname)
    plt.imshow(cutout_im_cutout, origin='lower')
    plt.savefig("/Users/edenmolina/Desktop/New PSF OutPut/%s_ChiSqCutout.pdf" % objectname)
    plt.close()

    Chi_Len = 25  # Number of possible normaliation constants to sample

    "FOR ITERATING THE VARIUS NORMALIZATION CONSTANTS"
    Multitplcation_Ratio = 5
    Norm_max = im1overim2 + im1overim2 * Multitplcation_Ratio
    Norm_min = im1overim2 - (im1overim2 * Multitplcation_Ratio / 10)
    x_Norm_space = np.linspace(Norm_min, Norm_max, Chi_Len)
    print "Org Ratio", im1overim2
    Chi_Sq_Array = np.zeros(len(x_Norm_space))

    "Explore Various Posible Normalizaition Constant"
    for i in range(len(x_Norm_space)):
        psf_subtracted = cutout_data - (cutout_im_cutout * x_Norm_space[i])
        Chi_Sq_Array[i] = np.sum(np.square(psf_subtracted))

    "Fitting the Data Aquired Above to a Quadratic, Then Minimize that Quadratic"
    x_quad_min, y_quad_min, constants = getMinimizedQuadratic(x_Norm_space, Chi_Sq_Array)

    "Plot the Chi Sq Plot"
    plotChiSq(x_Norm_space, constants, x_quad_min, y_quad_min, objectname, Chi_Sq_Array)
    plt.show()
    return x_quad_min

def saveImages(psf_sub, hdr, outdir, objectname):
    hdu = fits.PrimaryHDU(psf_sub, header=hdr)
    hdulist = fits.HDUList([hdu])
    hdulist.writeto("%s/%s_psfin.fits" % (outdir, objectname), clobber=True)

def getCutouts(image, obj1_coord, obj2_coord, simga_x):
    "Create cutouts of the images"
    # input coordinates as (y, x); in this case 16.58 is the sigma_x^2 of the object; radius is determined to be at Full Width Half Max
    #If Using circuilat citouts flip the ciirdinates using np.flip([], axis = 0)
    "Determines the size of the cutouts"
    image_separation = getZoom(obj1_coord, obj2_coord, image)
    cutout_im1 = square_mask(image, obj1_coord, 3.5*2 * np.sqrt(2 * np.log(2)) * simga_x[0])
    cutout_im2 = square_mask(image, obj2_coord,3.5*2 * np.sqrt(2 * np.log(2)) * simga_x[1])

    # Shifts for subtraction
    x1_shift, x2_shift, y1_shift, y2_shift, xobj_min, xobj_max, yobj_min, yobj_max = getshift(obj1_coord, obj2_coord)

    "Create the cutouts of the two lensed quasars and switch their positions; essentially shifts the images of the quasars"
    #This one places object 1 in the coordiantes of object 2
    im1_cutout = np.roll(
        np.roll(
            cutout_im1,
            int(x1_shift),
            axis=1),
        int(y1_shift),
        axis=0)

    im2_cutout = np.roll(
        np.roll(
            cutout_im2,
            int(x2_shift),
            axis=1),
        int(y2_shift),
        axis=0)

    #Returns the cutout of the two images, and then the shifted cutouts of the images
    return cutout_im1, cutout_im2, im1_cutout, im2_cutout

def cropCutout(im_cutout, obj1_coord, obj2_coord):
    "First find the min and max x y coordinates"
    x1_shift, x2_shift, y1_shift, y2_shift, xobj_min, xobj_max, yobj_min, yobj_max = getshift(obj1_coord, obj2_coord)

    "Get the maximum number of pixels from which we can extend away from one quasar coordinate in pixels"
    max_crop = getZoom(obj1_coord, obj2_coord, image)
    #Take a small fraction of this crop and use it to crop the images; Controls the zoom
    max_crop_frac = .2*max_crop

    xcrop_min = int(xobj_min - max_crop_frac)
    xcrop_max = int(xobj_max + max_crop_frac)
    ycrop_min = int(yobj_min - max_crop_frac)
    ycrop_max = int(yobj_max + max_crop_frac)
    im_cropped = im_cutout[ycrop_min:ycrop_max, xcrop_min:xcrop_max]

    # plt.close()
    # plt.subplot(121)
    # plt.title("Before")
    # plt.imshow(im_cutout, origin = 'lower')
    # plt.subplot(122)
    # plt.title("After")
    # plt.imshow(im_cropped, origin='lower')
    # plt.show()
    #Return the cropped image and the bounds of the image
    #saveImages(im_cropped, hdr, "/Users/edenmolina/Desktop/New PSF OutPut/", "S1128")

    return im_cropped, xcrop_min, ycrop_min

def getCroppedCoordiantes(obj_coord, xcrop_min, ycrop_min):
    "Returns the coordinates of the object based on the crop"
    x_new_crop = obj_coord[1] - xcrop_min
    y_new_crop = obj_coord[0] - ycrop_min

    new_coordiantes = np.array([y_new_crop, x_new_crop])
    return new_coordiantes

def plotShift(image, FileName, obj1_coord, obj2_coord, im1_cutout, im2_cutout, cutout_im1, cutout_im2):
    cmap=sns.cubehelix_palette(start=0.5, rot=-1.5, gamma=1., hue=1.,
                          # light=.8,
                          dark=0.,
                          reverse=True, as_cmap=True)
    from mpl_toolkits.axes_grid1 import make_axes_locatable #For the colorbar


    plt.close()
    "Plot the original data, and the two cutouts"
    fig=plt.figure(1, figsize=(18, 12))
    fig.subplots_adjust(hspace=.8)
    ax1 = plt.subplot(231)
    plt.title("%s Original Image" % FileName)
    image[image < 0.] = 1e-10
    im = plt.imshow(np.log10(image),  # use log10 as the scale than let's you compare change in magnitude
                     origin='lower',
                     vmin=-.9,
                     vmax=2,
                     cmap=sns.cubehelix_palette(start=0.5, rot=-1.5, gamma=1., hue=1.,
                                                # l/ight=.8,
                                                dark=0.,
                                                reverse=True, as_cmap=True))
    divider = make_axes_locatable(ax1)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = plt.colorbar(im, cax=cax)
    cbar.ax.set_ylabel(r'$\log_{10}({\rm counts/sec})$ (Flux)', rotation=90)

    ax2 = plt.subplot(232)
    plt.title("Object 1 Preshifted")
    plt.annotate('(%s, %s)' % ((obj1_coord[1], obj1_coord[0])),
                 xy=(obj1_coord[1], obj1_coord[0]), xytext=(image.shape[0] - 500, image.shape[1] - 120))
    im2 = plt.imshow(np.log10(cutout_im1),  # use log10 as the scale than let's you compare change in magnitude
                    origin='lower',
                    vmin=-.9,
                    vmax=2,
                    cmap=sns.cubehelix_palette(start=0.5, rot=-1.5, gamma=1., hue=1.,
                                               # l/ight=.8,
                                               dark=0.,
                                               reverse=True,  as_cmap=True))
    divider = make_axes_locatable(ax2)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = plt.colorbar(im2, cax=cax)
    cbar.ax.set_ylabel(r'$\log_{10}({\rm counts/sec})$ (Flux)', rotation=90)

    ax3=plt.subplot(233)
    plt.title("Object 2 Preshifted")
    plt.annotate('(%s, %s)' % ((obj2_coord[1], obj2_coord[0])),
                 xy=(obj2_coord[1], obj2_coord[0]), xytext=(image.shape[0] - 500, image.shape[1] - 120))

    im3 = plt.imshow(np.log10(cutout_im2),  # use log10 as the scale than let's you compare change in magnitude
                    origin='lower',
                    vmin=-.9,
                    vmax=2,
                    cmap=sns.cubehelix_palette(start=0.5, rot=-1.5, gamma=1., hue=1.,
                                               # l/ight=.8,
                                               dark=0.,
                                               reverse=True,  as_cmap=True))
    divider = make_axes_locatable(ax3)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = plt.colorbar(im3, cax=cax)
    cbar.ax.set_ylabel(r'$\log_{10}({\rm counts/sec})$ (Flux)', rotation=90)

    "Plot the shifted cutouts"
    ax4=plt.subplot(234)
    plt.title("%s Original Image" % FileName)
    image[image < 0.] = 1e-10
    im4 = plt.imshow(np.log10(image),  # use log10 as the scale than let's you compare change in magnitude
                    origin='lower',
                    vmin=-.9,
                    vmax=2,
                    cmap=sns.cubehelix_palette(start=0.5, rot=-1.5, gamma=1., hue=1.,
                                               # l/ight=.8,
                                               dark=0.,
                                               reverse=True,  as_cmap=True))
    divider = make_axes_locatable(ax4)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = plt.colorbar(im4, cax=cax)
    cbar.ax.set_ylabel(r'$\log_{10}({\rm counts/sec})$ (Flux)', rotation=90)

    ax5=plt.subplot(235)
    plt.title("Object 1 Shifted")
    im5 = plt.imshow(np.log10(im1_cutout),  # use log10 as the scale than let's you compare change in magnitude
                    origin='lower',
                    vmin=-.9,
                    vmax=2,
                    cmap=sns.cubehelix_palette(start=0.5, rot=-1.5, gamma=1., hue=1.,
                                               # l/ight=.8,
                                               dark=0.,
                                               reverse=True,  as_cmap=True))
    divider = make_axes_locatable(ax5)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = plt.colorbar(im5, cax=cax)
    cbar.ax.set_ylabel(r'$\log_{10}({\rm counts/sec})$ (Flux)', rotation=90)

    ax6=plt.subplot(236)
    plt.title("Object 2 Shifted")
    im6 = plt.imshow(np.log10(im2_cutout),  # use log10 as the scale than let's you compare change in magnitude
                     origin='lower',
                     vmin=-.9,
                     vmax=2,
                     cmap=sns.cubehelix_palette(start=0.5, rot=-1.5, gamma=1., hue=1.,
                                                # l/ight=.8,
                                                dark=0.,
                                                reverse=True, as_cmap=True))
    divider = make_axes_locatable(ax6)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = plt.colorbar(im6, cax=cax)
    cbar.ax.set_ylabel(r'$\log_{10}({\rm counts/sec})$ (Flux)', rotation=90)
    #plt.show()
    plt.savefig("/Users/edenmolina/Desktop/New PSF OutPut/%s_Shift.pdf" % FileName)
    plt.close()

def getInterpolation(image, res):
    "Interpolate the image"
    f = scinterp.interp2d(xrange(len(image[0])), xrange(len(image)), image)

    xx = np.arange(0, np.shape(image)[0], res)
    yy = np.arange(0, np.shape(image)[1], res)

    Interpolated_image = f(yy,xx)

    return Interpolated_image

#Invert the x-y coordiates
def getPSF(FileName, image, obj1_coord_org, obj2_coord_org, obj1_flux, obj2_flux, sigma_xsq, outputdir):
    "Set parameters"
    # Input parameters as [y, x]
    obj1_coord_org = obj1_coord_org
    obj2_coord_org = obj2_coord_org
    obj1_flux = obj1_flux
    obj2_flux = obj2_flux
    sigma_xsq = sigma_xsq
    flux1_over_flux2 = float(obj1_flux)/ float(obj2_flux)

    cutout_im1, cutout_im2, im1_cutout, im2_cutout = getCutouts(image, obj1_coord_org, obj2_coord_org, sigma_xsq)

    plt.close()
    plotShift(image, FileName, obj1_coord_org, obj2_coord_org, im1_cutout, im2_cutout, cutout_im1, cutout_im2)
    plt.close()
    "Get Chi Sq for the Two Objects"
    # Get getting Chi Square for Object 1, input the image cutout for the other object
    Norm_Obj1 = getChi_Sq(image, im2_cutout, flux1_over_flux2, np.flip(obj1_coord_org, axis=0),
                          np.array([0, 0]), FileName + 'Obj1')
    Norm_Obj2 = getChi_Sq(image, im1_cutout, flux1_over_flux2, np.flip(obj2_coord_org, axis=0),
                          np.array([0, 0]), FileName + 'Obj2')

    "CROP THE FILES"
    im1_cropped, xcrop_min, ycrop_min = cropCutout(im1_cutout, obj1_coord_org, obj2_coord_org)
    im2_cropped, xcrop_min, ycrop_min = cropCutout(im2_cutout, obj1_coord_org, obj2_coord_org)
    im_cropped, xcrop_min, ycrop_min = cropCutout(image, obj1_coord_org, obj2_coord_org)
    # Get new object coordinates
    obj1_coord = getCroppedCoordiantes(obj1_coord_org, xcrop_min, ycrop_min)
    obj2_coord = getCroppedCoordiantes(obj2_coord_org, xcrop_min, ycrop_min)

    InterpolationBool = False
    res = .1
    
    arcsecperpixel = .009
    v_min = -.9
    v_max = 1.8
    cmap=sns.cubehelix_palette(start=0.5, rot=-1.5, gamma=1., hue=1.,
                          # light=.8,
                          dark=0.,
                          reverse=True, as_cmap=True)
    plot_color = 'w'
    if InterpolationBool:

        plt.close()

        # "Interpolate the two images"
        print "Interpolating the two images"
        im1_inter = getInterpolation(im1_cropped, res)
        im2_im1_inter = getInterpolation(im2_cropped, res)
        im_inter = getInterpolation(im_cropped, res)
        print "Interpolation Complete"
        "Perform the PSF Subtraction"
        psf_subtracted = im_inter - im1_inter * Norm_Obj2 - im2_im1_inter * (Norm_Obj1)

        "Convert Pixels to Arcseconds"
        xsize_pixels = len(psf_subtracted[0])
        ysize_pixels = len(psf_subtracted)
        xsize_arcsec = xsize_pixels * arcsecperpixel * res
        ysize_arcsec = ysize_pixels * arcsecperpixel * res

        "Plot the PSF Subtraction"
        plt.figure(2, figsize=(14, 12))
        ax1 = plt.subplot(121)
        im = plt.imshow(np.log10(im_inter), origin='lower', cmap=cmap, extent=(-xsize_arcsec, xsize_arcsec, -ysize_arcsec, ysize_arcsec), vmin=v_min, vmax=v_max)
        # plt.plot(obj1_coord[1], obj1_coord[0], marker = 'o', color = 'g')
        # plt.plot(obj2_coord[1], obj2_coord[0], marker = 'o', color = 'g')
        plt.xlim(-xsize_arcsec, xsize_arcsec)
        plt.ylim(-ysize_arcsec, ysize_arcsec)
        plt.xlabel("arcsec")
        plt.ylabel("arcsec")
        plt.title("%s Observed" % FileName)
        #plt.colorbar()
        "Add the scale"
        xlim, ylim = ax1.get_xlim(), ax1.get_ylim()
        xstart = xlim[0] + 0.1 * (xlim[1] - xlim[0])  # This sets the start of the line 10% from left side of panel
        xend = xstart + 1  # This makes the line 1 arcsecond long
        y = ylim[0] + 0.05 * (ylim[1] - ylim[0])  # This sets the line 80% to the top of the panel
        ax1.plot([xstart, xend], [y, y], 'w-', lw=4, solid_capstyle='butt')  # Actually draw the line
        ax1.text(xstart + .5, y + .2, r'${1 \prime\prime}$', horizontalalignment='center', color='w')
        "North Direciton Arrow"
        x_arrow = xlim[0] + .9 * (xlim[1] - xlim[0])
        plt.arrow(x_arrow, y, 0, .5, head_width=.2, head_length=.3, fc='w', ec='w', width=.05)
        ax1.text(x_arrow, y + .9, 'N', horizontalalignment='center', color='w')
        plt.arrow(x_arrow, y, -.5, 0, head_width=.2, head_length=.3, fc='w', ec='w', width=.05)
        ax1.text(x_arrow - 1, y, 'E', horizontalalignment='center', color='w')


        ax2 = plt.subplot(122)
        plt.title("PSF Subtracted")
        plt.imshow(np.log10(psf_subtracted), origin='lower', cmap=cmap, vmin=v_min, vmax=v_max, extent=(-xsize_arcsec, xsize_arcsec, -ysize_arcsec, ysize_arcsec))
        # plt.plot(obj1_coord[1], obj1_coord[0], marker='o', color='g')
        # plt.plot(obj2_coord[1], obj2_coord[0], marker = 'o', color = 'g')
        #plt.colorbar()
        # Coordinates of the two objects
        plt.xlim(-xsize_arcsec, xsize_arcsec)
        plt.ylim(-ysize_arcsec, ysize_arcsec)
        plt.xticks([])
        plt.yticks([])
        "Plot the colorbar"
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        divider = make_axes_locatable(ax2)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        cbar = plt.colorbar(im, cax=cax)
        cbar.ax.set_ylabel(r'$\log_{10}({\rm (Flux)counts/sec})$', rotation=90)

        # plt.savefig("%s/%s_PSFSubtracted_Interpolated.pdf" % (outputdir, FileName))
        plt.savefig("%s/%s_PSFSubtracted_Interpolation.pdf" % (outputdir, FileName))
        print "PSF Completed", "\n"

    else:
        plt.close()
        "Perform the PSF Subtraction"
        psf_subtracted = im_cropped - im1_cropped * Norm_Obj2 - im2_cropped * (Norm_Obj1)

        "Convert Pixels to Arcseconds"
        # Coordinates of the two objects
        xsize_pixels = len(psf_subtracted[0])
        ysize_pixels = len(psf_subtracted)
        xsize_arcsec = xsize_pixels * arcsecperpixel
        ysize_arcsec = ysize_pixels * arcsecperpixel

        "Plot the PSF Subtraction"
        fig=plt.figure(2, figsize=(14, 12))
        fig.subplots_adjust(hspace=.8)
        ax1 = plt.subplot(121)
        im_cropped[im_cropped < 0.] = 1e-10
        im = plt.imshow(np.log10(im_cropped), origin='lower', cmap = cmap, extent=(-xsize_arcsec, xsize_arcsec, -ysize_arcsec, ysize_arcsec), vmin=v_min, vmax=v_max)
        # plt.plot(obj1_coord[1], obj1_coord[0], marker = 'o', color = 'g')
        # plt.plot(obj2_coord[1], obj2_coord[0], marker = 'o', color = 'g')
        plt.title("%s Observed" %FileName)
        #plt.colorbar()
        plt.xlim(-xsize_arcsec, xsize_arcsec)
        plt.ylim(-ysize_arcsec, ysize_arcsec)
        plt.xticks([])
        plt.yticks([])
        "Add the scale"
        xlim, ylim = ax1.get_xlim(), ax1.get_ylim()
        xstart = xlim[0] + 0.1 * (xlim[1] - xlim[0])  # This sets the start of the line 10% from left side of panel
        xend = xstart + 1  # This makes the line 1 arcsecond long
        y = ylim[0] + 0.05 * (ylim[1] - ylim[0])  # This sets the line 80% to the top of the panel
        ax1.plot([xstart, xend], [y, y], '%s-'%plot_color, lw=4, solid_capstyle='butt')  # Actually draw the line
        ax1.text(xstart + .5, y + .2, r'${1 \prime\prime}$', horizontalalignment='center', color=plot_color)
        "North Direciton Arrow"
        x_arrow = xlim[0] + .9 * (xlim[1] - xlim[0])
        plt.arrow(x_arrow, y, 0, .5, head_width=.2, head_length=.3, fc=plot_color, ec=plot_color, width=.05)
        ax1.text(x_arrow, y + .9, 'N', horizontalalignment='center', color=plot_color)
        plt.arrow(x_arrow, y, -.5, 0, head_width=.2, head_length=.3, fc=plot_color, ec=plot_color, width=.05)
        ax1.text(x_arrow - 1, y, 'E', horizontalalignment='center', color=plot_color)

        ax2 = plt.subplot(122)
        plt.title("PSF Subtracted")
        psf_subtracted[psf_subtracted < 0.] = 1e-10
        im=plt.imshow(np.log10(psf_subtracted), origin='lower', cmap = cmap, extent=(-xsize_arcsec, xsize_arcsec, -ysize_arcsec, ysize_arcsec), vmin=v_min, vmax=v_max)
        # plt.plot(obj1_coord[1], obj1_coord[0], marker='o', color='g')
        # plt.plot(obj2_coord[1], obj2_coord[0], marker = 'o', color = 'g')
        plt.xticks([])
        plt.yticks([])
        "Plot the colorbar"
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        divider = make_axes_locatable(ax2)
        #cax = divider.append_axes("right", size="5%", pad=.05)
        plt.subplots_adjust(bottom=0.1, right=0.8, top=0.8)
        cax = plt.axes([0.85, 0.15, 0.030, 0.6])
        cbar = plt.colorbar(im, cax=cax)
        cbar.ax.set_ylabel(r'$\log_{10}({\rm (Flux)counts/sec})$', rotation=90)


        #plt.savefig("%s/%s_PSFSubtracted_Interpolated.pdf" % (outputdir, FileName))
        plt.savefig("%s/%s_PSFSubtracted_NoInterpolation.pdf" %(outputdir,FileName))

        print "PSF Completed"
        print "Saved to ", "%s/%s_PSFSubtracted_NoInterpolation.pdf" %(outputdir,FileName)

        return psf_subtracted

def saveImages(psf_sub,hdr,outdir,objectname):
    hdu = fits.PrimaryHDU(psf_sub, header=hdr)
    hdulist = fits.HDUList([hdu])
    hdulist.writeto("%s%s_psfsubtracted.fits" % (outdir,objectname),clobber=True)
    print "%s%s_psfsubtracted.fits" % (outdir,objectname)

if __name__ == "__main__":
    "Import the image data"
    #Data = 'October'
    #Data = 'January'
    #Data = 'Rerun'
    Data = '2013'
    #Data = '2016'
    #Data = 'Rerun'
    PSF_input = pd.read_excel("/Users/edenmolina/Documents/Astro/PSFSubtraction/PSF_Input.xlsx", sheetname=Data)
    FileNames = PSF_input['Object Name']
    print "Found ", len(FileNames), " Objects"
    X1 = PSF_input['x1']
    Y1 = PSF_input['y1']
    Fluxes1 = PSF_input['f1']
    X2 = PSF_input['x2']
    Y2 = PSF_input['y2']
    Fluxes2 = PSF_input['f2']
    SigmaX1 = PSF_input['sigx1']
    SigmaX2 = PSF_input['sigx2']
    ImagePath = PSF_input['filepath']

    for i in range(len(FileNames)):
        #hdulist = fits.open("/Users/edenmolina/Documents/Astro/Coadds/%s2018Coadds/%s_coadd.fits" %(Data,FileNames[i]))
        hdulist = fits.open(ImagePath[i])
        hdr = hdulist[0].header
        image = hdulist[0].data

        image -= np.median(image)  # Remove background

        print "Working on ", FileNames[i]

        "Set parameters"
        #Input parameters as [y, x]
        obj1_coord_org = np.array([Y1[i], X1[i]])
        obj2_coord_org = np.array([Y2[i], X2[i]])
        obj1_flux = Fluxes1[i]
        obj2_flux = Fluxes2[i]
        sigma_xsq = [SigmaX1[i], SigmaX2[i]]

        psf = getPSF(FileNames[i], image, obj1_coord_org, obj2_coord_org, obj1_flux, obj2_flux, sigma_xsq, "/Users/edenmolina/Desktop/New PSF OutPut")
        saveImages(psf, hdr, "/Users/edenmolina/Desktop/New PSF OutPut/", FileNames[i])
