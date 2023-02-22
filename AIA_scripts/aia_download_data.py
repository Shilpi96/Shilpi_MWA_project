import astropy.units as u
from sunpy.net import Fido, attrs
import sunpy.map
import os
import glob
#from aiapy.calibrate import register, update_pointing, normalize_exposure
# 193 angstorm --- 1 & 20 millikelvin - slightly hotter and very hot material of a solar flare
# 171 angstorm --- 600000 kelvin -- upper transition region / quiet corona
# 211 angstorm --- 2 million kelvin -- active regions

wavelngth = 193
result = Fido.search(attrs.Time('2014/08/25 15:04:30', '2014/08/25 15:05'),
                     attrs.Instrument.aia,
                     attrs.Wavelength(wavelngth*u.angstrom))
print(result)
files = Fido.fetch(result, path = '/mnt/LOFAR-PSP/STELLAR-20140825-EVENT/AIA_data/aiaimages_'+str(wavelngth)+'/', overwrite = False)
print('saved the fit files')

list1 = glob.glob('/mnt/LOFAR-PSP/STELLAR-20140825-EVENT/AIA_data/aiaimages_'+str(wavelngth)+'/*.fits')

###### Now we have to delete files for which the exposure time is less than the standard time. for 171 = 1.5 sec, 211 = 1.5 sec and 193 = 1.0 sec. because during solar eruptions the image can get saturated so the camera will try to close more quickly because it will try to receive less photons. That means we are not getting correct information. So need to delete these images.

for i in range(len(list1)):
	aia = sunpy.map.Map(list1[i])
	if aia.meta['exptime']<1.0 : os.remove(list1[i])

print('deleted the fit files which had the exposure time less than 1.5 sec')
