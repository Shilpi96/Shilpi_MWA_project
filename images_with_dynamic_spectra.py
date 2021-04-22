import matplotlib.pyplot as plt
import numpy as np

from astropy import units as u
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.io import fits
from astropy.time import Time
import pdb
import sunpy.data.sample
import sunpy.map
from sunpy.coordinates import frames, sun
from scipy import optimize
from astropy.modeling.functional_models import Gaussian2D
from astropy.io import fits
from astropy.wcs import WCS
from pylab import *
from numpy import *    ### it is needed to define ravel
import datetime
from datetime import timedelta
from matplotlib.ticker import MaxNLocator
from matplotlib import dates



##############################################################################
## Get the visibility data
a = np.load('/mnt/LOFAR-PSP/shilpi_MWA_project/1095907576/dynamic_spectrum_from_interferometric_data/visibilities.npy')
b = np.load('/mnt/LOFAR-PSP/shilpi_MWA_project/CLEAN1_1095907576/1095907576_103-104/fitsfile/dynamic_spectrum_from_interferometric_data/visibilities.npy')
a = np.delete(a, slice(586,592), axis = 3)
b = np.delete(b, slice(582,592), axis = 3)

vis  = np.concatenate((a,b), axis = 3)
  
freq = np.load('/mnt/LOFAR-PSP/shilpi_MWA_project/1095907576/dynamic_spectrum_from_interferometric_data/freq.npy')


print('Calculating times ...')
start_time = datetime.datetime(2014, 9, 28, 2, 46, 6, 0)
timestep = timedelta(seconds = 0.5)
times = [(start_time + i*timestep) for i in range(len(vis[0][0][0]))]
times_mpl = [dates.date2num(t) for t in times]
'''ylabels = [[131.5, 132.5,133.5], [118.5, 119.5, 120.5], 
            [107.0, 108.0, 109.0], [97.0, 98.0, 99.0], 
            [88.0, 89.0, 90.0], [79.0, 80.0, 81.0]]'''
ylabels = [[131.5, 133.5], [118.5,  120.5], 
            [107.0,  109.0], [97.0,  99.0], 
            [88.0,  90.0], [79.0,  81.0]]

######## Creating the submap
list1 = ['/mnt/LOFAR-PSP/shilpi_MWA_project/CLEAN1_1095907576/1095907576_103-104/fitsfile/103-104/Sun_Clean_Matrix_025000.0_53~55_full_fc.fits','/mnt/LOFAR-PSP/shilpi_MWA_project/CLEAN1_1095907576/1095907576_103-104/fitsfile/093-094/Sun_Clean_Matrix_025000.0_53~55_full_fc.fits','/mnt/LOFAR-PSP/shilpi_MWA_project/CLEAN1_1095907576/1095907576_103-104/fitsfile/084-085/Sun_Clean_Matrix_025000.0_53~55_full_fc.fits','/mnt/LOFAR-PSP/shilpi_MWA_project/CLEAN1_1095907576/1095907576_103-104/fitsfile/076-077/Sun_Clean_Matrix_025000.0_53~55_full_fc.fits','/mnt/LOFAR-PSP/shilpi_MWA_project/CLEAN1_1095907576/1095907576_103-104/fitsfile/069-070/Sun_Clean_Matrix_025000.0_53~55_full_fc.fits','/mnt/LOFAR-PSP/shilpi_MWA_project/CLEAN1_1095907576/1095907576_103-104/fitsfile/062-063/Sun_Clean_Matrix_025000.0_53~55_full__fc.fits']
Lofar_submap = []
Frequency = []
for i in range(len(list1)):
	hdu = fits.open(list1[i])[0]
	header = hdu.header
	data = hdu.data[0, 0, :, :]
	obstime = Time(header['date-obs'])
	frequency = header['crval3']*u.Hz
	Frequency.append(frequency)
	#pdb.set_trace()
	lofar_loc = EarthLocation(lat=-26.7033*u.deg, lon=116.671*u.deg)
	lofar_gcrs = SkyCoord(lofar_loc.get_gcrs(obstime))
	reference_coord = SkyCoord(header['crval1']*u.Unit(header['cunit1']),
                           header['crval2']*u.Unit(header['cunit2']),
                           frame='gcrs',
                           obstime=obstime,
                           obsgeoloc=lofar_gcrs.cartesian,
                           obsgeovel=lofar_gcrs.velocity.to_cartesian(),
                           distance=lofar_gcrs.hcrs.distance)
	reference_coord_arcsec = reference_coord.transform_to(frames.Helioprojective(observer=lofar_gcrs))
	cdelt1 = (np.abs(header['cdelt1'])*u.deg).to(u.arcsec)
	cdelt2 = (np.abs(header['cdelt2'])*u.deg).to(u.arcsec)
	P1 = sun.P(obstime)
	new_header = sunpy.map.make_fitswcs_header(data, reference_coord_arcsec,
                                           reference_pixel=u.Quantity([header['crpix1']-1,
                                                                       header['crpix2']-1]*u.pixel),
                                           scale=u.Quantity([cdelt1, cdelt2]*u.arcsec/u.pix),
                                           rotation_angle=-P1,
                                           wavelength=frequency.to(u.MHz).round(2),
                                           observatory='LOFAR')
	lofar_map = sunpy.map.Map(data, new_header)
	lofar_map_rotate = lofar_map.rotate()
	bl = SkyCoord(-3000*u.arcsec, -2500*u.arcsec, frame=lofar_map_rotate.coordinate_frame)
	tr = SkyCoord(3000*u.arcsec, 2500*u.arcsec, frame=lofar_map_rotate.coordinate_frame)
	lofar_submap = lofar_map_rotate.submap(bl, top_right=tr)
	Lofar_submap.append(lofar_submap)
###### Making the plots
fig = plt.figure(figsize=(12, 9))
grid = plt.GridSpec(21, 3, wspace=0.2, hspace=0.25)
ax1 = plt.subplot(grid[0, :3])     ###### the plot will span 1st to 3rd column
ax2 = plt.subplot(grid[1, :3])
ax3 = plt.subplot(grid[2, :3])
ax4 = plt.subplot(grid[3, :3])
ax5 = plt.subplot(grid[4, :3])
ax6 = plt.subplot(grid[5, :3])
ax7 = plt.subplot(grid[6:13, 0],projection = Lofar_submap[0]) ##### the plot will span row 6th to 13th

Lofar_submap[0].plot(cmap='jet')
Lofar_submap[0].draw_limb(color='darkmagenta')
Lofar_submap[0].draw_grid(color='gray')
ax7.set_xlabel('  ')
ax7.set_ylabel('  ')
ax7.set_title('  ')
ax7.tick_params(axis='x',labelbottom = False)
text(0.95, 0.05, """Freq(MHz) : %.1f """        #### (0.95,0.05) is the position of the text on the plot
%(Frequency[0].value/10**6),fontsize=7, horizontalalignment='right',
verticalalignment='bottom', transform=ax7.transAxes, color = 'white')
text(0.95, 0.03,obstime.value[11:]+" UT" , horizontalalignment='right',
				verticalalignment='bottom', transform=ax7.transAxes,fontsize = 7, color ="white")

ax8 = plt.subplot(grid[6:13, 1],projection = Lofar_submap[1])
Lofar_submap[1].plot(cmap='jet')
Lofar_submap[1].draw_limb(color='darkmagenta')
Lofar_submap[1].draw_grid(color='gray')
ax8.set_xlabel('  ')
ax8.set_ylabel('  ')
ax8.tick_params(axis='x',labelbottom = False)
ax8.tick_params(axis='y',labelbottom = False)
ax8.set_title('  ')

ax9 = plt.subplot(grid[6:13, 2],projection = Lofar_submap[2])
Lofar_submap[2].plot(cmap='jet')
Lofar_submap[2].draw_limb(color='darkmagenta')
Lofar_submap[2].draw_grid(color='gray')
ax9.set_xlabel('  ')
ax9.set_ylabel('  ')
ax9.tick_params(axis='x',labelbottom = False)
ax9.tick_params(axis='y',labelbottom = False)
ax9.set_title('  ')

ax10 = plt.subplot(grid[14:21, 0],projection = Lofar_submap[3])
Lofar_submap[3].plot(cmap='jet')
Lofar_submap[3].draw_limb(color='darkmagenta')
Lofar_submap[3].draw_grid(color='gray')
ax10.set_xlabel('  ')
ax10.set_ylabel('  ')
ax10.set_title('  ')

ax11 = plt.subplot(grid[14:21, 1],projection = Lofar_submap[4])
Lofar_submap[4].plot(cmap='jet')
Lofar_submap[4].draw_limb(color='darkmagenta')
Lofar_submap[4].draw_grid(color='gray')
ax11.set_xlabel('  ')
ax11.set_ylabel('  ')
ax11.tick_params(axis='y',labelbottom = False)
ax11.set_title('  ')

ax12 = plt.subplot(grid[14:21, 2],projection = Lofar_submap[5])
Lofar_submap[5].plot(cmap='jet')
Lofar_submap[5].draw_limb(color='darkmagenta')
Lofar_submap[5].draw_grid(color='gray')
ax12.set_xlabel('  ')
ax12.set_ylabel('  ')
ax12.tick_params(axis='y',labelbottom = False)
ax12.set_title('  ')
axs = [ax1,ax2,ax3,ax4,ax5,ax6]


########## Plotting the Lofar_submap in loop
'''for i in range(6):
	Lofar_submap[i].plot(cmap='jet')
	Lofar_submap[i].draw_limb(color='darkmagenta')
	Lofar_submap[i].draw_grid(color='gray')'''

######## plotting the dynamic spectra
vmin = np.log10(1000)
vmax = np.log10(180000)
t0 = dates.date2num(times[0])
t1 = dates.date2num(times[-1])

for i in range(6):
    XX = vis[i][0]*np.conj(vis[i][0])
    YY = vis[i][3]*np.conj(vis[i][3])
    f0 = freq[i][0][0]
    f1 = freq[i][-1][0]
    spectra = np.log10(np.sqrt(np.real(XX + YY))/2)
    peak = axs[i].imshow(spectra, 
    cmap=plt.get_cmap('Spectral_r'), 
    vmin=vmin, vmax = vmax, aspect = 'auto', origin = 'lower',
    extent=(t0, t1, f0, f1))
    axs[i].set_yticks(ylabels[i])
    axs[i].set_yticklabels(ylabels[i])
    axs[i].tick_params(axis='y', labelsize= 9)
    axs[i].xaxis.tick_top()
    axs[i].axvline(x = times[468], color = 'black')
    #axs[i].label_outer()
ax1.xaxis.tick_top()
ax1.xaxis_date()
ax1.xaxis.set_major_formatter(dates.DateFormatter('%H:%M:%S'))
ax1.set_xlim(times_mpl[0], times_mpl[len(vis[0][0][0])-1])
#ax6.set_xticklabels(fontsize = 6)
ax6.tick_params(axis='x', labelsize= 9)
for i in range(1,6):
	axs[i].label_outer()
	axs[i].set_xticklabels(' ')
	axs[i].tick_params(axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off         # ticks along the top edge are off
    labelbottom=False)
plt.show()

