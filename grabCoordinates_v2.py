"""
Stable Version of grabCoordinates
To be used with MCMC_CircularGaussian

This is the GUI promt that will ask the user to choose an image file (.fits format) and specify which part of the image to calculate the flux of a bright object such as a quasar or galaxy
    Depends on MCMC_CircularGaussian to run. This contains the Monte Carlo Markov Chain parameters and algorithm
"""
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
import astropy.io.fits as fits
#import MCMCTest_v5 as MCMC
import MCMC_CircularGaussian as MCMC
from matplotlib.colors import LogNorm
import Tkinter as Tk
import tkFileDialog

class ImagePositions:
    def __init__(self):
        self.coords = []
        # This is to get things set up so that mouse clicks and other events
        # in matplotlib will be registered
        self.cid = fig.canvas.mpl_connect('button_press_event', self)

    def __call__(self, event):
        # When you click on the image, print off the coordinates that were
        # clicked and append to coords.
        print "Image at (%.2f , %.2f)" % (event.xdata,event.ydata)
        self.coords.append([event.xdata,event.ydata])

"Select the File with a GUI Promt"
root = Tk.Tk()
root.withdraw()
file_path = tkFileDialog.askopenfilename()
file_save_name = raw_input("Please Input The Object Name: ")
#file_save_name = file_path[-21:-11]
print "Working on: ", file_save_name

"File path of test"
root.destroy()

hdulist = fits.open(file_path)
data = hdulist[0].data

"Plot the image data; Displays in Log Scale"
fig, ax = plt.subplots(figsize = (10, 10))
zoom_level = 0 #Set how much to zoom in for selecting the images
ax.imshow(data[int(np.min(data[0])+zoom_level):int(np.max(data[0])-zoom_level),
          int(np.min(data[1])+zoom_level):int(np.max(data[1])-zoom_level)],
          cmap="gray", norm = LogNorm(), vmin = 1, vmax = 25)
ax.imshow(data, cmap="gray", norm = LogNorm(), vmin = 1, vmax = 25, origin = 'lower')

imagepositions = ImagePositions()

"Display the selected file"
plt.show()
plt.close()
"Print Selected Coordinates to User"
Coordinates = np.array(imagepositions.coords)

"Run MCMC For each of the selected coordinates"
sample = list()

"Running the MCMC Code for each one of the cooridnates selected"
for i in range(len(Coordinates)):
    x =  Coordinates[i][0]
    y = Coordinates[i][1]
    "Grabs the flux at given coordinate"
    Flux_Guess = data[int(round(y))][int(round(x))]

    zoom = 50 #Set how many pizels to process around the selected coordinate
    #plt.imshow(data[int(y - zoom):int(y + zoom), int(x - zoom):int(x + zoom)], origin = 'lower')
    #plt.imshow(data[int(y - zoom):int(y + zoom), int(x - zoom):int(x + zoom)], origin='lower', cmap = "gray")
    #plt.show()
    "Crops the image data to be within the 'zoom' area; speeds up processing"
    sampledata = MCMC.Flux(file_save_name, i ,data[int(y - zoom):int(y + zoom), int(x - zoom):int(x + zoom)], Flux_Guess, zoom, zoom, (x-zoom), (y-zoom))
