import numpy as np
import pdb
import glob
import sunpy.map
import pylab
from datetime import datetime
from astropy.io import fits
from astropy.time import Time
from astropy import units as u
import matplotlib.pyplot as plt
from scipy import stats
import matplotlib.dates as mdates
import matplotlib
matplotlib.rcParams.update({'font.size': 11}) 

### Get the distance 

def Dis(gauss_file, fint=0, tint=0):
	pos = gauss_file	

	xasec = pos[fint][tint][0]
	yasec = pos[fint][tint][1]
	
	
	Coords = [354, -348] #### Position of the active region
	dis = np.sqrt(np.square(xasec - Coords[0]) + np.square(yasec - Coords[1]))
	b = 696340/957.57 ### 1 arcsec = 767 km
	dis = dis*b/1e3
	return dis
	
### Function to get the index number of the points

def get_index(radio_times, lrmnth_times, i = 2):
	delt = abs(radio_times - lrmnth_times[:, np.newaxis])
	index = []
	for i in range(i):
		delt[i] = np.where(delt[i] == delt[i].min())
		index.append(delt[i][0][0])
	return index
	

##### Load the time and position of the points from AIA
pos_20 = np.load('/mnt/LOFAR-PSP/shilpi_MWA_project/AIA_data/npy_files/pos_20_EUV_front.npy')
pos_40 = np.load('/mnt/LOFAR-PSP/shilpi_MWA_project/AIA_data/npy_files/pos_40_EUV_front.npy')
pos_60 = np.load('/mnt/LOFAR-PSP/shilpi_MWA_project/AIA_data/npy_files/pos_60_EUV_front.npy')
pos_80 = np.load('/mnt/LOFAR-PSP/shilpi_MWA_project/AIA_data/npy_files/pos_80_EUV_front.npy')
pos_100 = np.load('/mnt/LOFAR-PSP/shilpi_MWA_project/AIA_data/npy_files/pos_100_EUV_front.npy')
pos_120 = np.load('/mnt/LOFAR-PSP/shilpi_MWA_project/AIA_data/npy_files/pos_120_EUV_front.npy')
pos_140 = np.load('/mnt/LOFAR-PSP/shilpi_MWA_project/AIA_data/npy_files/pos_140_EUV_front.npy')

time = np.load('/mnt/LOFAR-PSP/shilpi_MWA_project/AIA_data/npy_files/Time_EUV_front.npy')

#pdb.set_trace()
##### Calculate the distance for each of the point from the active region
b = 696340/957.57   ### to calculate how much is 1 arcsecond in km on the Sun. radius = 957.57 arcsec = 696340 km
Coords = [354, -348]  ### pos of the active region
Pos_20 = np.zeros(38) 
Pos_40 = np.zeros(38) 
Pos_60 = np.zeros(38) 
Pos_80 = np.zeros(38) 
Pos_100 = np.zeros(38) 
Pos_120 = np.zeros(38) 
Pos_140 = np.zeros(38)
 
for i in range(38):
	
	a = np.sqrt(np.square(pos_20[i][0] - Coords[0]) + np.square(pos_20[i][1] - Coords[1]))
	Pos_20[i] = a*b/1e3  #### 1 arsecond = 767 km and converting it to megameter

	a = np.sqrt(np.square(pos_40[i][0] - Coords[0]) + np.square(pos_40[i][1] - Coords[1]))
	Pos_40[i] = a*b/1e3
	
	a = np.sqrt(np.square(pos_60[i][0] - Coords[0]) + np.square(pos_60[i][1] - Coords[1]))
	Pos_60[i] = a*b/1e3
	
	a = np.sqrt(np.square(pos_80[i][0] - Coords[0]) + np.square(pos_80[i][1] - Coords[1]))
	Pos_80[i] = a*b/1e3
	
	a = np.sqrt(np.square(pos_100[i][0] - Coords[0]) + np.square(pos_100[i][1] - Coords[1]))
	Pos_100[i] = a*b/1e3
	
	a = np.sqrt(np.square(pos_120[i][0] - Coords[0]) + np.square(pos_120[i][1] - Coords[1]))
	Pos_120[i] = a*b/1e3
	
	a = np.sqrt(np.square(pos_140[i][0] - Coords[0]) + np.square(pos_140[i][1] - Coords[1]))
	Pos_140[i] = a*b/1e3

##### Getting the time
time = Time(time)
aia_time = [datetime.strptime(Time(time)[i].value[:-1], "%Y-%m-%dT%H:%M:%S.%f") for i in range(len(time))]	

aia_time = [datetime.timestamp(aia_time[i]) for i in range(len(aia_time))]
aia_time  = np.array(aia_time)
AIA_time = aia_time-aia_time[0]

##### Performing linear regression
res_20 = stats.linregress(AIA_time, Pos_20)
res_40 = stats.linregress(AIA_time, Pos_40)
res_60 = stats.linregress(AIA_time, Pos_60)
res_80 = stats.linregress(AIA_time, Pos_80)
res_100 = stats.linregress(AIA_time, Pos_100)
res_120 = stats.linregress(AIA_time, Pos_120)
res_140 = stats.linregress(AIA_time, Pos_140)

#pdb.set_trace()
##### Plotting the line and points
fig = plt.figure(figsize=(8, 8))

## colorscale for plotting the pointa and the line other that the 60 degrees line
cm = pylab.get_cmap('Blues')
collist = cm(np.linspace(0, 255, 25).astype(int))

plt.plot(AIA_time,Pos_20,'o',color = collist[7])
plt.plot(AIA_time, res_20.intercept + res_20.slope*AIA_time, color = collist[7], label = '20 deg, {:.2f}'.format((res_20.slope)*1e3)+' km/s')
plt.plot(AIA_time,Pos_40,'o', color = collist[14])
plt.plot(AIA_time, res_40.intercept + res_40.slope*AIA_time, color = collist[12], label = '40 deg, {:.2f}'.format((res_40.slope)*1e3)+' km/s')
plt.plot(AIA_time,Pos_60,'o', color = 'Green')
plt.plot(AIA_time, res_60.intercept + res_60.slope*AIA_time, color = 'Green', label = '60 deg, {:.2f}'.format((res_60.slope)*1e3)+' km/s')
plt.plot(AIA_time,Pos_80,'o', color = collist[19])
plt.plot(AIA_time, res_80.intercept + res_80.slope*AIA_time, color = collist[17], label = '80 deg, {:.2f}'.format((res_80.slope)*1e3)+' km/s')

#plt.plot(AIA_time,Pos_100,'o',color = 'red')
#plt.plot(AIA_time, res_100.intercept + res_100.slope*AIA_time, 'red', label = '100 deg, {:.2f}'.format((res_100.slope)*1e3)+' km/s')
plt.plot(AIA_time,Pos_120,'o',color = collist[20])
plt.plot(AIA_time, res_120.intercept + res_120.slope*AIA_time, color = collist[24], label = '120 deg, {:.2f}'.format((res_120.slope)*1e3)+' km/s')
#plt.plot(AIA_time,Pos_140,'o',color = 'purple')
#plt.plot(AIA_time, res_140.intercept + res_140.slope*AIA_time, 'purple', label = '140 deg, {:.2f}'.format((res_140.slope)*1e3)+' km/s')'''

#### Plotting LFB radio sources
root = '/mnt/LOFAR-PSP/shilpi_MWA_project/'

x_point_LFB = np.load(root+'AIA_data/npy_files/x_coord_LFB.npy')
x_point_HFB = np.load(root+'AIA_data/npy_files/x_coord_HFB.npy')

### How to convert datetime objects to datetime.datetime objects

time = mdates.date2num(x_point_LFB) ### convert it to matplotlib dates objects
time = mdates.num2date(x_point_LFB) ### then convert it to datetime.datetime objects

LFB_time = [time[0],time[7], time[8], time[14], time[15],time[22],time[23],time[30],time[31],time[38]]
LFB_time = [LFB_time[i].replace(tzinfo=None) for i in range(len(LFB_time))]
LFB_time = np.array(LFB_time)

time1 = mdates.date2num(x_point_HFB) ### convert it to matplotlib dates objects
time1 = mdates.num2date(x_point_HFB) ### then convert it to datetime.datetime objects

HFB_time = [time1[0],time1[7], time1[8], time[15], time[16],time[23],time[27],time[30]]

HFB_time = [HFB_time[i].replace(tzinfo=None) for i in range(len(HFB_time))]
HFB_time = np.array(HFB_time)

## Load the time from the radio files for both the datasets
time_7576_103 = np.load(root+'AIA_data/npy_files/7576_radio_time/7576_103-104_radio_time.npy', allow_pickle=True)
time_7576_093 = np.load(root+'AIA_data/npy_files/7576_radio_time/7576_093-094_radio_time.npy', allow_pickle=True)
time_7576_084 = np.load(root+'AIA_data/npy_files/7576_radio_time/7576_084-085_radio_time.npy', allow_pickle=True)
time_7576_076 = np.load(root+'AIA_data/npy_files/7576_radio_time/7576_076-077_radio_time.npy', allow_pickle=True)
time_7872_084 = np.load(root+'AIA_data/npy_files/7872_radio_time/7872_084-085_radio_time.npy', allow_pickle=True)
time_7872_076 = np.load(root+'AIA_data/npy_files/7872_radio_time/7872_076-077_radio_time.npy', allow_pickle=True)
time_7872_069 = np.load(root+'AIA_data/npy_files/7872_radio_time/7872_069-070_radio_time.npy', allow_pickle=True)

### Get the index for chosen time for LFB
Lindex_103 = get_index(time_7576_103, LFB_time[0:2])
Lindex_093 = get_index(time_7576_093, LFB_time[2:4])
Lindex_084 = get_index(time_7576_084, LFB_time[4:6])
Lindex_076 = get_index(time_7576_076, LFB_time[6:8])
Lindex_069 = get_index(time_7872_069, LFB_time[8:10])

### Get the index for chosen time for HFB

Hindex_103 = get_index(time_7576_103, HFB_time[0:2])
Hindex_093 = get_index(time_7576_093, HFB_time[2:4])
Hindex_084 = get_index(time_7872_084, HFB_time[4:6])
Hindex_076 = get_index(time_7872_076, HFB_time[6:8])

### Upload the npy files containig the informations of the Gaussian fitted sources
gauss_file0 = np.load(root+'CLEAN1_1095907576/new_analysis/103-104/coords_in_arcsec.npy')
gauss_file1 = np.load(root+'CLEAN1_1095907576/new_analysis/093-094/coords_in_arcsec.npy')
Gauss_file0 = np.load(root+'CLEAN1_1095907576/new_analysis/084-085/coords_in_arcsec.npy')
Gauss_file1 = np.load(root+'CLEAN1_1095907576/new_analysis/076-077/coords_in_arcsec.npy')
Gauss_file2 = np.load(root+'CLEAN1-1095907872/new_analysis/069-070/coords_in_arcsec.npy')
Gauss_file3 = np.load(root+'CLEAN1-1095907872/new_analysis/084-085/coords_in_arcsec.npy')
Gauss_file4 = np.load(root+'CLEAN1-1095907872/new_analysis/076-077/coords_in_arcsec.npy')


### Get the distance for the source positions using the index for LFB
Ldis = []
Ldis[0:2] = [Dis(gauss_file0, fint=5, tint=Lindex_103[0]),Dis(gauss_file0, fint=7, tint=Lindex_103[1]) ]
Ldis[2:4] = [Dis(gauss_file1, fint=0, tint=Lindex_093[0]),Dis(gauss_file1, fint=7, tint=Lindex_093[1]) ]
Ldis[4:6] = [Dis(Gauss_file0, fint=0, tint=Lindex_084[0]),Dis(Gauss_file0, fint=7, tint=Lindex_084[1]) ]
Ldis[6:8] = [Dis(Gauss_file1, fint=0, tint=Lindex_076[0]),Dis(Gauss_file1, fint=7, tint=Lindex_076[1]) ]
Ldis[8:10] = [Dis(Gauss_file2, fint=0, tint=Lindex_069[0]),Dis(Gauss_file2, fint=7, tint=Lindex_069[1]) ]

Ldis = np.array(Ldis)

### Get the distance for the source positions using the index for HFB
Hdis = []
Hdis[0:2] = [Dis(gauss_file0, fint=0, tint=Hindex_103[0]),Dis(gauss_file0, fint=7, tint=Hindex_103[1]) ]
Hdis[2:4] = [Dis(gauss_file1, fint=0, tint=Hindex_093[0]),Dis(gauss_file1, fint=7, tint=Hindex_093[1]) ]
Hdis[4:6] = [Dis(Gauss_file3, fint=0, tint=Hindex_084[0]),Dis(Gauss_file3, fint=7, tint=Hindex_084[1]) ]
Hdis[6:8] = [Dis(Gauss_file4, fint=0, tint=Hindex_076[0]),Dis(Gauss_file4, fint=7, tint=Hindex_076[1]) ]

Hdis = np.array(Hdis)

### Get the time for the sources in second and do linear regression
LFB_Time = np.zeros(10)
for i in range(10):
	LFB_Time[i] = datetime.timestamp(LFB_time[i]) 
radio_time = LFB_Time-aia_time[0]
res1 = stats.linregress(radio_time, Ldis)

### Get the time for the sources in second and do linear regression for HFB
HFB_Time = np.zeros(HFB_time.shape[0])
for i in range(HFB_time.shape[0]):
	HFB_Time[i] = datetime.timestamp(HFB_time[i]) 
radio_Time = HFB_Time-aia_time[0]
res2 = stats.linregress(radio_Time, Hdis)


plt.plot(radio_time,Ldis,'v',color = 'green',markersize = 8, mec = 'black',label = 'LFB, {:.2f}'.format((res1.slope)*1e3)+' km/s')
plt.plot(radio_time, res1.intercept + res1.slope*radio_time, 'green')
'''
plt.plot(radio_Time,Hdis,'o',color = 'blue')
plt.plot(radio_Time, res2.intercept + res2.slope*radio_Time, 'blue', label = 'HFB, {:.2f}'.format((res2.slope)*1e3)+' km/s')'''

plt.ylabel('Distance (Mm)')
plt.xlabel('Time (sec)')
plt.legend(loc = 'upper right')
#plt.xlim(140,500)
plt.show()
#pdb.set_trace()
