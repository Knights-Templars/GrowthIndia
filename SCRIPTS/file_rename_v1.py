import os
import glob
import shutil
from astropy.io import fits

filelist=sorted(glob.glob('*.fits'))
def file_rename(files):
	for i,filename in enumerate(files,1):
		hdul=fits.open(filename)
		OBJECT=hdul[0].header['OBJECT'].replace(" ","")[0:6]
		FILTER=hdul[0].header['FILTER'].replace(" ","")[-1]
		EXPTIME=hdul[0].header['EXPTIME']
                if OBJECT=='Biasim':
                    new_name=str(OBJECT)+"_"+str(i)+".fits"
                else:
                    new_name=str(OBJECT)+"_"+str(FILTER)+"_"+str(i)+".fits"
		os.rename(filename,new_name)
file_rename(filelist)


