import os
import glob
import pdb
#########################################################################
# Part of the scripts needed to align all of the 20140928 TypeII images
# for Shilpi's project.
#
# This script does the following:
# - Use CASA tclean to make shallow cleaned images for each of the fine
#   channels all of the coarse channels and MS
#
# TO DO FOR SHILPI
# - Confirm if IMSIZE, CELL are what were used for the earlier imaging.
#
#   Divya 22Jun2021
#########################################################################

#Loop on MS
#Loop on coarse channels in a MS

def define_TIMES(MSNAME,TIME_STEP):
#returns a list of time slices, taking steps defined in TIME_STEPS 
    ms.open(MSNAME)
    ms.select({'antenna1':[1],'antenna2':[1]})
    tt = ms.getdata('time')
    ms.close()
    times = tt['time'][::TIME_STEP]
    dt = float(tt['time'][1])-float(tt['time'][0])
    
    print ('#### First time stamp %s' % (tt['time'][0]))
    print ('#### Integration time %f' % (dt))
    print ('#### Number of time stamps %d' % (len(tt['time'])))
    print ('#### Last time stamp %s' % (tt['time'][len(tt['time'])-1]))
    print ('#### No. of time stamps to image %d' % (len(times)))
    return times, dt

MS_DIR='/mnt/LOFAR-PSP/shilpi_MWA_project'
GAINTABLE_DIR='/mnt/LOFAR-PSP/shilpi_MWA_project'
#MS_LIST = ['CLEAN1_1095907576', 'CLEAN1-109507872']
MS_LIST = ['CLEAN1_1095907576']
#COARSE_CHAN_LIST = ['062-063', '069-070', '076-077', '084-085', '093-094', '103-104']
COARSE_CHAN_LIST = ['093-094'] 
#FINE_CHAN_LIST = ['chan_8~11', 'chan_12~15', 'chan_17~20', 'chan_21~23', 'chan_40~43', 'chan_44~47', 'chan_49~52', 'chan_53~55']
FINE_CHAN_LIST = ['chan_12~15']

IMSIZE = 1620
CELL = '100arcsec'
NITER = 200

TIME_STEP = 1

for MS in MS_LIST:
	OBS_ID = '1095'+(MS.split('1095')[1])
	print (OBS_ID)

	for COARSE_CHAN in COARSE_CHAN_LIST:
		j = 0
		for FINE_CHAN in FINE_CHAN_LIST:
			MS1 = MS_DIR+'/'+MS+'/MS/'+COARSE_CHAN+'/'+FINE_CHAN+'.ms'
			if os.path.isdir(MS1) == 0: # Check if MS is available
				print (MS1)

#			Figure out the timerange, but needs to be done only once for each COARSE_CHANNEL
			if (j == 0):
				TIME_LIST, DT = define_TIMES(MS1,TIME_STEP)
			IMAGE_DIR = '/mnt/LOFAR-PSP/shilpi_MWA_project/'+MS+'/New_Images/093_094'
			IMAGE_NAME_BASE = IMAGE_DIR+'/'+OBS_ID+'_'+COARSE_CHAN+'_'+FINE_CHAN
			i = 0
			while (i < len(TIME_LIST)):
#			while (i < 101):
				index = '%04d' % (i*TIME_STEP)

				TIME0=qa.time(qa.quantity(TIME_LIST[i]-DT/2.0,'s'),form="ymd")
				TIME1=qa.time(qa.quantity(TIME_LIST[i]+DT/2.0,'s'),form="ymd")
				time0=str(TIME0).split('\'')[1]+'.0' # The split is get rid of the '''s
				time1=str(TIME0).split('\'')[1]+'.5'
				TIME_IN=time0+'~'+time1
				temp = time0.split('/')[3].split(':')
				TIME_TAG = temp[0]+temp[1]+temp[2]
				IMAGENAME=IMAGE_NAME_BASE+'_'+TIME_TAG
				print (IMAGENAME)

				print ("### Cleaning MS %s Time = %s" %(MS1, TIME_IN))
				tclean(vis=MS1,selectdata=True,field="",spw="",
					timerange=TIME_IN,uvrange="",antenna="",scan="",observation="",intent="",datacolumn="corrected",
					imagename=IMAGENAME,
					imsize=IMSIZE,
					cell=CELL,phasecenter="",stokes="I",projection="SIN",startmodel="",specmode="mfs",reffreq="",nchan=-1,start="",width="",outframe="LSRK",	veltype="radio",restfreq=[],interpolation="linear",perchanweightdensity=False,gridder="standard",facets=1,psfphasecenter='',chanchunks=1,wprojplanes=1,vptable="",usepointing=False,mosweight=True,aterm=True,psterm=False,wbawp=True,conjbeams=False,cfcache="",computepastep=360.0,rotatepastep=360.0,pblimit=0.2,normtype="flatnoise",deconvolver="hogbom",scales=[],nterms=1,smallscalebias=0.6,restoration=True,restoringbeam=[],pbcor=False,outlierfile="",weighting="natural",robust=0.5,noise="1.0Jy",npixels=0,uvtaper=[],niter=NITER,gain=0.1,threshold=0.0,nsigma=0.0,cycleniter=-1,cyclefactor=1.0,minpsffraction=0.05,maxpsffraction=0.8,interactive=False,usemask="user",mask="",pbmask=0.0,sidelobethreshold=3.0,noisethreshold=5.0,lownoisethreshold=1.5,negativethreshold=0.0,smoothfactor=1.0,minbeamfrac=0.3,cutthreshold=0.01,growiterations=75,dogrowprune=True,minpercentchange=-1.0,verbose=False,fastnoise=True,restart=True,savemodel="none",calcres=True,calcpsf=True,parallel=False)
				i = i + 1
			j = j + 1
