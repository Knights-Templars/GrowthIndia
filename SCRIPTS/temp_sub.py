import ois
import matplotlib.pyplot as plt
from astropy.io import fits
ref_image='reference_image_005.resamp.fits'
sci_image='20190428142418-005-RA.wcs.resamp.fits'
ref_image_open=fits.open(ref_image)
sci_image_open=fits.open(sci_image)
ref_data=ref_image_open[0].data
sci_data=sci_image_open[0].data
#diff=ois.optimal_system(sci_data, ref_data,method='Alard-Lupton',gausslist=[{'sx': 0.1, 'sy':0.1, 'modPolyDeg':0}],bkgdegree=2,kernelshape=(11,11))[0]
diff=ois.optimal_system(sci_data, ref_data,method='AdaptiveBramich',poly_degree=2)[0]
plt.imshow(diff,cmap='viridis',interpolation='none',origin='lower')
hdu_sub_image=fits.PrimaryHDU(diff)
hdu_sub_image.writeto('sub_final_11.fits')
