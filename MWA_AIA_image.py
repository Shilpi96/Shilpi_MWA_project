from sunpy.coordinates.sun import sky_position as sun_position
import sunpy.coordinates.sun as sun_coord
import numpy as np
import datetime
from datetime import timedelta
import sunpy.map
import matplotlib.pyplot as plt
import astropy.units as unit
import pdb
import pylab
from astropy.coordinates import SkyCoord
import matplotlib.colors as colors
from astropy.visualization import ImageNormalize,LogStretch

####### loading the RA data and converting it from radian to degree
pos1 = np.load('/Users/shilpibhunia/Documents/MWA_Project/Data/Helioprojective_coordinates/RA_DEC_103.npy')
pos2 = np.load('/Users/shilpibhunia/Documents/MWA_Project/Data/Helioprojective_coordinates/RA_DEC_093.npy')
pos3 = np.load('/Users/shilpibhunia/Documents/MWA_Project/Data/Helioprojective_coordinates/RA_DEC_084.npy')
pos4 = np.load('/Users/shilpibhunia/Documents/MWA_Project/Data/Helioprojective_coordinates/RA_DEC_076.npy')

##### Download and load the AIA data

map1 = sunpy.map.Map('/Users/shilpibhunia/Documents/MWA_Project/Data/aia_02:46:13_171_lev1.5.fits')
top_right = SkyCoord(1250 * unit.arcsec, 400 * unit.arcsec, frame=map1.coordinate_frame)
bottom_left = SkyCoord(-350 * unit.arcsec, -1350 * unit.arcsec, frame=map1.coordinate_frame)
Map1 = map1.submap(bottom_left, top_right=top_right)
fig = plt.figure()
ax = fig.add_subplot(111, projection=Map1)
Map1.plot(cmap = 'magma',norm=ImageNormalize(vmin=0,vmax=5000,stretch=LogStretch(10)))
cm = pylab.get_cmap('YlOrRd')
cm1 = pylab.get_cmap('Greens')
cm2 = pylab.get_cmap('Blues')
cm3 = pylab.get_cmap('cool')
cm4 = pylab.get_cmap('copper')
cm5 = pylab.get_cmap('BuPu')
pdb.set_trace()
X = []
Y = []
for i in range(443):
    X.append(pos1[7][i][0]*unit.arcsec.to('deg'))
    Y.append(pos1[7][i][1]*unit.arcsec.to('deg'))
   
#pdb.set_trace()
'''for i in range(443):
    color = cm(i)
    ax.plot(X[i:i+1],Y[i:i+1],linestyle='-')
num_lines = 443

colors = [plt.cm.Blues(i) for i in np.linspace(0, 1, num_lines)]

ax.plot(X,Y, 'o',transform=ax.get_transform('world'), color = colors)'''

#ax.plot(X, Y,'-',transform=ax.get_transform('world'), color = 'red')
for i in range(443):
    color = cm(i)
    ax.plot(pos1[0][i][0]* unit.arcsec.to('deg'), pos1[0][i][1] * unit.arcsec.to('deg'), 'o', color = color, transform=ax.get_transform('world'))
ax.plot(X,Y,'-',transform=ax.get_transform('world'),color = 'red')
'''for i in range(ntimes1):
    color1 = cm1(i)
    ax.plot(y[i][0]* unit.arcsec.to('deg'), y[i][1] * unit.arcsec.to('deg'), 'o', color = color1, transform=ax.get_transform('world'),
        label='120.56 MHz')
for i in range(ntimes2):
    color2 = cm2(i)
    ax.plot(z[i][0]* unit.arcsec.to('deg'), z[i][1] * unit.arcsec.to('deg'), 'o', color = color2, transform=ax.get_transform('world'),
        label='109.04 MHz')
for i in range(ntimes3):
    color3 = cm3(i)
    ax.plot(l[i][0]* unit.arcsec.to('deg'), l[i][1] * unit.arcsec.to('deg'), 'o', color = color3, transform=ax.get_transform('world'),
        label='98.80 MHz')
for i in range(ntimes4):
    color4 = cm4(i)
    ax.plot(k[i][0]* unit.arcsec.to('deg'), k[i][1] * unit.arcsec.to('deg'), 'o', color = color4, transform=ax.get_transform('world'),
        label='89.84 MHz')
for i in range(ntimes5):
    color5 = cm5(i)
    ax.plot(uu[i][0]* unit.arcsec.to('deg'), uu[i][1] * unit.arcsec.to('deg'), 'o', color = color5, transform=ax.get_transform('world'),
        label='80.88 MHZ')'''

Map1.draw_grid()
plt.show()