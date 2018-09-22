"""
Stable Version

A file dialog will pop up and promt the user to choose an image file (.fits format)
A new window displaying the image will show and user will select area on image where user wants to learn about the flux at that spot
Takes a defined window size given by variable size on Line 53 and computes that root mean square, mean, and total flux of a selection.

"""
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
import astropy.io.fits as fits
from matplotlib.colors import LogNorm
import Tkinter as Tk
import tkFileDialog
import pandas as pd

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
file_save_name = file_path[-21:-11]
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
#print Coordinates[0][0], Coordinates[0][1]


size = 50

x_10 = np.linspace(int(Coordinates[0][0])-(size/2), int(Coordinates[0][0])+((size/2)), size)
y_10 = np.linspace(int(Coordinates[0][1])-(size/2), int(Coordinates[0][1])+((size/2)), size)

"Run MCMC For each of the selected coordinates"
FluxArray = np.zeros((size, size))

"Running the MCMC Code for each one of the cooridnates selected"
for i in range(len(x_10)):
    x = x_10[i]
    for j in range(len(y_10)):
        y = y_10[j]
        "Grabs the flux at given coordinate"
        Image_Value = data[int(y)][int(x)]
        FluxArray[i][j] = Image_Value


Avg = np.average(FluxArray)
print "The Root Mean Square of the 100: ", np.sqrt(np.sum(np.mean(FluxArray**2), axis = None))
print "The Mean Flux of the 100 is: ", Avg
print "The Total Count is: ", np.sum(FluxArray)
print "With Average Subtracted: ", np.sum(FluxArray-Avg)


plt.close()
figure = plt.figure(figsize=(14, 10))
plot1 = plt.subplot(121)
plot1.hist(data[int(y_10[0]): int(y_10[-1]), int(x_10[0]): int(x_10[-1])].flatten(), bins=19)
plt.xlabel("Pixel Value")
plt.ylabel("Frequency")
#plt.hist(data[int(y_10[0]): int(y_10[-1]), int(x_10[0]): int(x_10[-1])], bins=20)
print np.shape(data[int(y_10[0]): int(y_10[-1]), int(x_10[0]): int(x_10[-1])])
plot2 = plt.subplot(122)
plot2.imshow(data[int(y_10[0]): int(y_10[-1]), int(x_10[0]): int(x_10[-1])])
plt.xlabel("X")
plt.ylabel("Y")
plt.show()