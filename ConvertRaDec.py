"""
Converts right assention and declanation from HH:MM:SS to Degrees
Will be used to determine the names of the observed systems
"""

import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u

def ConvertRaDecFormatted(ra='', dec='', round=False):
    RA, DEC, rs, ds = '', '', '', ''

    if dec:
        if str(dec)[0] == '-':
            ds, dec = '-', abs(dec)
        deg = int(dec)
        decM = abs(int((dec - deg) * 60))
        if round:
            decS = int((abs((dec - deg) * 60) - decM) * 60)
        else:
            decS = (abs((dec - deg) * 60) - decM) * 60

        if ds == '-':
            DEC = '-'+ str(deg).zfill(2) + str(decM).zfill(2)
        else:
            DEC = '+' + str(deg).zfill(2) + str(decM).zfill(2)
        #DEC = str(ds).zfill(2) + str(deg).zfill(2)

    if ra:
        if str(ra)[0] == '-':
            rs, ra = '-', abs(ra)
        raH = int(ra / 15)
        raM = int(((ra / 15) - raH) * 60)
        if round:
            raS = int(((((ra / 15) - raH) * 60) - raM) * 60)
        else:
            raS = ((((ra / 15) - raH) * 60) - raM) * 60
        if rs == '-':
            RA = str(raH).zfill(2) + str(raM).zfill(2)
        else:
            RA = str(raH).zfill(2) + str(raM).zfill(2)

    if ra and dec:
        return ("J"+RA+DEC)
    else:
        return ("J" + RA+ DEC)



RADEC = np.loadtxt("/Users/edenmolina/PycharmProjects/Quasar/RaDec/RADEC2016.txt")

for i, ra in enumerate(RADEC[0]):
    "Convert using Astropy"
    #Astropy_Result = SkyCoord(ra=ra*u.degree, dec=RADEC[1][i]*u.degree)


    "Print the Conversions"
    #print "Using own code: ", ConvertRaDecFormatted(ra=ra, dec=RADEC[1][i], round=False)
    try:
        print ConvertRaDecFormatted(ra=ra, dec=RADEC[1][i], round=False)
    except:
        print ra, RADEC[1][i]
    #print "Using Astropy: ", "J", int(Astropy_Result.ra.hms[0]), int(Astropy_Result.ra.hms[1]), Astropy_Result.dec.degree