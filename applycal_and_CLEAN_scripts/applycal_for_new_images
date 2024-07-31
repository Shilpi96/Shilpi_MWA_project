
import os
import glob

#########################################################################
# Part of the scripts needed to align all of the 20140928 TypeII images
# for Shilpi's project.
#
# This script does the following:
# - Applycal solutions corresponding to a single timestamp to all of the data (2 MS)
# - Needs to be done for each of the fine channels, for each of the coarse channels
#
# - It turns out that the solutions are not available for all of the time stamps
# - So one has to choose the times (CAL_TIME_LIST) for each fine channel, coarse
#   channel and MS.
# 
# - To avoid overwriting the original MS, copy over all of the MS corresponding for
#   each of the fine and coarse channels for both the MS to an independent area
#   before applycal.
#
#  TO DO FOR SHILIPI before running this script
# - Copy over data for CLEAN1_1095907576 for COARSE CHANNELS '076-077', '084-085', '093-094' and '103-104'
#   to MS_DIR in the same format as it is present there.
# - Copy over data for CLEAN1-1095907872 for all COARSE CHANNELS to MS_DIR
# - Make sure to provide a list of time stamps for which solutions are available 
#   in the gaintable directory (e.g. /data1/shilpi/CLEAN1_1095907576/1095907576_103-104/chan_40~43)
#
#   Divya 22Jun2021
#########################################################################

#Loop on MS
#Loop on coarse channels in a MS

MS_DIR='/mnt/LOFAR-PSP/shilpi_MWA_project'
GAINTABLE_DIR='/mnt/LOFAR-PSP/shilpi_MWA_project'
#MS_LIST = ['CLEAN1_1095907576', 'CLEAN1-1095907872']
MS_LIST = ['CLEAN1_1095907576']
#COARSE_CHAN_LIST = ['062-063', '069-070', '076-077', '084-085', '093-094', '103-104']
COARSE_CHAN_LIST = ['093-094'] # A shorter list to try out when testing the code
#FINE_CHAN_LIST = ['chan_8~11', 'chan_12~15', 'chan_17~20', 'chan_21~23', 'chan_40~43', 'chan_44~47', 'chan_49~52', 'chan_53~55']
FINE_CHAN_LIST = ['chan_12~15'] # A shorter list to try out when testing the code

for MS in MS_LIST:
	OBS_ID = '1095'+(MS.split('1095')[1])
	print OBS_ID

	for COARSE_CHAN in COARSE_CHAN_LIST:
		i = 0
		if (OBS_ID == '1095907576'):
			'''if (COARSE_CHAN == '062-063'):
				CAL_TIME_LIST = ['024900.0', '024900.0', '024900.0', '024900.0', '024900.0', '024900.0', '024900.0', '024900.0']
			if (COARSE_CHAN == '069-070'):
				CAL_TIME_LIST = ['024648.0', '024648.0', '024648.0', '024648.0', '024648.0', '024648.0', '024648.0', '024648.0']
			if (COARSE_CHAN == '076-077'):
				CAL_TIME_LIST = ['024648.0', '024648.0', '024648.0', '024648.0', '024648.0', '024648.0', '024648.0', '024648.0']
			if (COARSE_CHAN == '084-085'):
				CAL_TIME_LIST = ['024648.0', '024648.0', '024648.0', '024648.0', '024648.0', '024648.0', '024648.0', '024648.0']
			if (COARSE_CHAN == '093-094'):
				CAL_TIME_LIST = ['024606.0', '024606.0', '024606.0', '024606.0', '024606.0', '024606.0', '024606.0', '024606.0']
			if (COARSE_CHAN == '103-104'):
				CAL_TIME_LIST = ['025009.0', '025009.0', '025009.0', '025009.0', '025009.5', '025009.0', '025010.0', '025009.0']'''
			

		'''if (OBS_ID == '1095907872'):
			if (COARSE_CHAN == '062-063'):
				CAL_TIME_LIST = ['025116.0', '025116.0', '025116.5', '025116.0', '025116.5', '025114.0', '025120.0', '025116.0']
			if (COARSE_CHAN == '069-070'):
				CAL_TIME_LIST = ['025537.5', '025537.5', '025524.5', '025537.5', '025537.5', '025537.5', '025543.5', '025537.5']
			if (COARSE_CHAN == '076-077'):
				CAL_TIME_LIST = ['025548.0', '025548.0', '025548.0', '025548.0', '025548.0', '025539.5', '025548.0', '025547.5']
			if (COARSE_CHAN == '084-085'):
				CAL_TIME_LIST = ['025506.0', '025506.0', '025506.0', '025506.0', '025506.0', '025506.0', '025506.0', '025506.0']
			if (COARSE_CHAN == '093-094'):
				CAL_TIME_LIST = ['025506.0', '025506.0', '025506.0', '025506.0', '025506.0', '025506.0', '025506.0', '025506.0']

			if (COARSE_CHAN == '103-104'):
				CAL_TIME_LIST = ['025102.0', '025506.0', '025506.5', '025506.0', '025506.5', '025506.5', '025506.0', '025506.5']'''


		for FINE_CHAN in FINE_CHAN_LIST:
			MS1 = MS_DIR+'/'+MS+'/MS/'+COARSE_CHAN+'/'+FINE_CHAN+'.ms'
			print(MS1)
			GAINTABLE = GAINTABLE_DIR+'/'+MS+'/'+OBS_ID+'_'+COARSE_CHAN+'/'+FINE_CHAN+'/caltables/caltable_'+CAL_TIME_LIST[i]+'.cal'
			print(GAINTABLE)
			if os.path.isdir(MS1) == 0:
				print MS1
			if os.path.isdir(GAINTABLE) == 0:
				print GAINTABLE
#applycal
			print "### Running clearcal on MS %s" %(MS1)
			clearcal(vis=MS1,field="",spw="",intent="",addmodel=False)
			print "### Applying gaintable %s to MS %s" %(GAINTABLE, MS1)
			applycal(vis=MS1,field="",spw="",intent="",selectdata=True,timerange="",uvrange="",antenna="",scan="",observation="",msselect="",docallib=False,callib="",gaintable=GAINTABLE,gainfield=[],interp=[],spwmap=[],calwt=[True],parang=False,applymode="",flagbackup=False)

			i = i + 1


