import glob
import os
import numpy as np
from astropy.io import fits


filelist = glob.glob('*.fits')

for f, file in enumerate(filelist):
	img = fits.open(filelist[f])
	ra = img[0].header['TARRA']
	dec = img[0].header['TARDEC']
	img.close()
	#print(ra)
	#print(dec)
	print(str(filelist[f]))
	os.system('solve-field'+" "+" --ra "+str(ra)+","+" --dec "+str(dec)+","+" --radius 1.0"+","+" --crpix-center"+" "+"--tweak-order 3"+","+" --overwrite"+" "+" --sigma 5"+","+" --nsigma 20.0"+","+" --no-background-subtraction"+" --resort"+" "+" --keep-xylist %s.xy"+" "+str(filelist[f]))
	print("The given command is:")
	#print(astrometry_cmd)
	#status=os.system(astrometry_cmd)
	#print(status)
newlist = glob.glob('*.new')
print('Following files have been generated. Remaining files will be deleted'%(newlist))

for file in newlist:
    file_wcs=file.replace("new", "wcs.fits")
    print(file_wcs)
    os.system("mv "+file+" "+file_wcs)



