# A Python Code for performing Alignment of Astronomical Images having WCS 
# Author - Anirban Dutta, Avinash Singh
# Version - 1.0
# Date: 04/07/2019

#---------------------------------------------------------------------------------------------------------------------#
# Import the necessary modules and packages

import os
import numpy as np
import shutil
import glob
from astropy.io import fits
import subprocess
import numpy
from astropy.stats import sigma_clipped_stats
from astropy.wcs import WCS
import warnings
from astroquery.gaia import Gaia
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.io import ascii


from pyraf import iraf

#---------------------------------------------------------------------------------------------------------------------#

# Test whether SCAMP, SWARP and SEXTRACTOR are installed properly

dependencies=['swarp', 'scamp', 'sextractor']

def test_dependencies(dep):
    try:
       subprocess.Popen(dep, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
       print("%s is installed properly. OK" % dep)
       return 1
    except ImportError:
        print("==%s IS NOT INSTALLED PROPERLY" % dep)
        return 0
        
i=0
for dep in dependencies:
    i+=test_dependencies(dep)
print("\n%i out of %i dependencies installed properly." % (i, len(dependencies)))
if i!=len(dependencies):
    print('**********Please Install the programs before continuing**********')
else:
    print('**********You are ready to continue**********')
    
#---------------------------------------------------------------------------------------------------------------------#

# Telescope and CCD Specifications:

Telescope_name= 'GROWTH_INDIA'
CCD_name='RTS2'
read_noise=14
gain=1.6
data_max=55000
OBJECT='OBJECT'
OBJECT_NAME='SN_2019ein'
Right_Ascension='RA'
Declination='DEC'
RA='13:53:29.134'
DEC='+40:16:31.4'

#---------------------------------------------------------------------------------------------------------------------#

cwd=os.getcwd()
config_dir=cwd+"/CONFIG/"
DIR_aligned=cwd+"/SN_ALIGNED/"

#---------------------------------------------------------------------------------------------------------------------#
# function for removing files

def remove_file(file_name):
    try:
        os.remove(file_name)
    except OSError:
        pass
        
# function for removing files having similar names

def remove_similar_files(common_text):
    for residual_file in glob.glob(common_text):
        remove_file(residual_file)

def group_similar_files(text_list, common_text, exceptions=''):
    list_files=glob.glob(common_text)
    if exceptions !='' :
        list_exceptions=exceptions.split(',')
        for file_name in list_files:
            for text in list_exceptions:
                test=re.search(text, file_name)
                if test:
                    list_files.remove(file_name)

    list_files.sort()
    if len(text_list) !=0:
        with open(text_list, 'w') as f:
            for file_name in list_files:
                f.write(file_name+'\n')
                
    return list_files    
    


def text_list_to_python_list(text_list):
    if os.path.exists(text_list):
        with open(text_list, 'r+') as f:
            python_list=f.read.split()
            
    return python_list
    
def python_list_to_text_list(python_list, text_list):
    with open(text_list, 'w') as f:
        for element in python_list:
            f.write(str(element)+'\n')
            
            
            
def list_lists_to_list(list_lists, text_list):
    list_name=[]
    for file_name in list_lists:
        with open(file_name, 'r') as f:
            file_list=f.read().split()
            for element in file_list:
                list_name.append(element)
    python_list_to_text_list(list_name, text_list)
    
    return list_name    
    
#---------------------------------------------------------------------------------------------------------------------#


for text in ['*.xyls', '*.axy', '*.corr', '*.match', '*.new', '*.wcs', '*.solved', '*.rdls', '*.png', 'wcs_*', 
         '*.list', "*.resamp.fits", "*.resamp.weight.fits", "log_*", "list_*", '*.ldac', '*.txt', 'log*']:

    remove_similar_files(common_text=text)

iraf.noao(_doprint=0)
iraf.images(_doprint=0)

# function to edit the header

def hedit(textlist_files, field_name, value, add_keyword='no'):

    task=iraf.images.imutil.hedit
    task.unlearn()
    
    task.add=str(add_keyword)
    task.addonly='no'
    task.delete='no'
    task.verify='no'
    task.show='no'
    task.update='yes'
    task.mode='ql'
    
    task(images='@'+textlist_files, fields=field_name, value=value)            
   
# -------------------------------------------------------------------------------------------------------#

# Ignoring warnings in output 

warnings.filterwarnings("ignore")

# -------------------------------------------------------------------------------------------------------#

"""
FUNCTIONS TO CONVERT FITS FILES OR ASTROPY TABLES TO FITS_LDAC FILES AND
VICE VERSA.

"""

def convert_hdu_to_ldac(hdu):
    
    """
    Convert an hdu table to a fits_ldac table
   
    """
    from astropy.io import fits
    import numpy as np
    tblhdr = np.array([hdu.header.tostring(',')])
    col1 = fits.Column(name='Field Header Card', array=tblhdr, format='13200A')
    cols = fits.ColDefs([col1])
    tbl1 = fits.BinTableHDU.from_columns(cols)
    tbl1.header['TDIM1'] = '(80, {0})'.format(len(hdu.header))
    tbl1.header['EXTNAME'] = 'LDAC_IMHEAD'
    tbl2 = fits.BinTableHDU(hdu.data)
    tbl2.header['EXTNAME'] = 'LDAC_OBJECTS'
    return (tbl1, tbl2)

def convert_table_to_ldac(tbl):
    
    """
    Convert an astropy table to a fits_ldac
    
    """
    from astropy.io import fits
    import tempfile
    f = tempfile.NamedTemporaryFile(suffix='.fits', mode='rb+')
    tbl.write(f, format='fits')
    f.seek(0)
    hdulist = fits.open(f, mode='update')
    tbl1, tbl2 = convert_hdu_to_ldac(hdulist[1])
    new_hdulist = [hdulist[0], tbl1, tbl2]
    new_hdulist = fits.HDUList(new_hdulist)
    return new_hdulist

def save_table_as_ldac(tbl, filename, **kwargs):
    
    """
    Save a table as a fits LDAC file
    
    """
    hdulist = convert_table_to_ldac(tbl)
    hdulist.writeto(filename, **kwargs)



def get_table_from_ldac(filename, frame=1):
    
    """
    Load an astropy table from a fits_ldac by frame 
    
    """
    from astropy.table import Table
    if frame>0:
        frame = frame*2
    tbl = Table.read(filename, hdu=frame)
    return tbl

# -------------------------------------------------------------------------------------------------------#

            
# function to run Astrometry.net and find initial estimate of WCS of the images
  
def run_astrometry(filelist):

    list_files=glob.glob(filelist)
    for filename in list_files:
        hdul=fits.open(filename)
        RA=hdul[0].header['TARRA']
        DEC=hdul[0].header['TARDEC']
        OBJECT=hdul[0].header['OBJECT']
        data_image=hdul[0].data
        mean, median, std_dev=sigma_clipped_stats(data_image, sigma=3.0)
        nsigma=20*median
        
        astrometry_command='solve-field'+" "+" --ra "+str(RA)+","+" --dec "+str(DEC)+"," \
        +" --radius 1.0"+","+" --crpix-center"+" "+" --tweak-order 2"+","+" --overwrite"+" "+"--resort" \
        +" --new-fits"+" "+"wcs_"+filename+" "+filename
        
        print("----------Astrometry.net-->Solve-Field is running----------")
        print(astrometry_command)
        os.system(astrometry_command)      
        
        
        
'''
Function to query gaia catalog with a given RA, DEC and a search radius

'''

def make_gaia_catalog(ra, dec, catalog_box_size, catalog_min_mag, catalog_max_mag, catname, writetext = True, writeldac = True):
	job = Gaia.launch_job_async("SELECT * FROM gaiadr2.gaia_source AS g, gaiadr2.panstarrs1_best_neighbour AS pbest, gaiadr2.panstarrs1_original_valid AS ps1 WHERE g.source_id = pbest.source_id AND pbest.original_ext_source_id = ps1.obj_id AND CONTAINS(POINT('ICRS', g.ra, g.dec), CIRCLE('ICRS', %.4f, %.4f, %.4f))=1 AND ps1.r_mean_psf_mag > %.2f AND ps1.r_mean_psf_mag < %.2f AND pmra IS NOT NULL AND pmdec IS NOT NULL AND abs(pmdec) > 0 AND abs(pmdec) < 40 AND abs(pmra)>0 AND abs(pmra) < 40 AND ps1.n_detections > 6 AND pbest.number_of_mates=0 AND pbest.number_of_neighbours=1;"%(ra, dec, catalog_box_size, catalog_min_mag, catalog_max_mag), dump_to_file = False)
	
	p = job.get_results()
	
	# convert RA and DEC errors from mas(milli arc second) to degrees

	p['ra_errdeg'] = p['ra_error'] / 3.6e6
	p['dec_errdeg'] = p['dec_error'] / 3.6e6
	p['FLAGS']=0
	
	
	p.remove_columns(['astrometric_n_obs_al', 'astrometric_n_obs_ac', 'astrometric_n_good_obs_al', 'astrometric_n_bad_obs_al', 'astrometric_gof_al', 'astrometric_chi2_al', 'astrometric_excess_noise', 'astrometric_excess_noise_sig', 'astrometric_params_solved', 'astrometric_primary_flag', 'astrometric_weight_al', 'astrometric_pseudo_colour', 'astrometric_pseudo_colour_error', 'mean_varpi_factor_al', 'astrometric_matched_observations', 'visibility_periods_used', 'astrometric_sigma5d_max', 'frame_rotator_object_type', 'matched_observations', 'duplicated_source', 'phot_g_n_obs', 'phot_g_mean_flux', 'phot_g_mean_flux_error', 'phot_g_mean_flux_over_error', 'phot_g_mean_mag', 'phot_bp_n_obs', 'phot_bp_mean_flux', 'phot_bp_mean_flux_error', 'phot_bp_mean_flux_over_error', 'phot_bp_mean_mag', 'phot_rp_n_obs', 'phot_rp_mean_flux', 'phot_rp_mean_flux_error', 'phot_rp_mean_flux_over_error', 'phot_rp_mean_mag', 'phot_bp_rp_excess_factor', 'phot_proc_mode', 'bp_rp', 'bp_g', 'g_rp', 'radial_velocity', 'radial_velocity_error', 'rv_nb_transits', 'rv_template_teff', 'rv_template_logg', 'rv_template_fe_h', 'l', 'b', 'ecl_lon', 'ecl_lat', 'priam_flags', 'teff_val', 'teff_percentile_lower', 'teff_percentile_upper', 'a_g_val', 'a_g_percentile_lower', 'a_g_percentile_upper', 'e_bp_min_rp_val', 'e_bp_min_rp_percentile_lower', 'e_bp_min_rp_percentile_upper', 'flame_flags', 'radius_val', 'radius_percentile_lower', 'radius_percentile_upper', 'lum_val', 'lum_percentile_lower', 'lum_percentile_upper', 'gaia_astrometric_params', 'obj_name', 'obj_id', 'ra_2', 'dec_2', 'ra_error_2', 'dec_error_2', 'epoch_mean', 'zone_id', 'obj_info_flag', 'quality_flag', 'designation', 'phot_variable_flag', 'datalink_url', 'epoch_photometry_url', 'original_ext_source_id'])
	
	if writetext:
	    if os.path.exists(catname+'.txt'):
	        os.remove(catname+'.txt')
	    ascii.write(p, catname+'.txt')
	if writeldac:
		if os.path.exists(catname + '.ldac'):
			os.remove(catname + '.ldac')
		save_table_as_ldac(p, catname + '.ldac')
    
    
'''
Run Sextractor on the list of images and create thier catalogs

'''
def run_sextractor(ctext):

    for text in ['*.list']:
        remove_similar_files(common_text=text)

    cwd=os.getcwd()
    config_sex=cwd+'/config.sex'
    param_sex=cwd+'/sextractor.param'

    file_list='ImageList.list'
    list_files=group_similar_files(file_list, common_text=ctext)

    
    
    
    for file_name in list_files:
        command="sextractor %s -c %s -CATALOG_NAME %s -PARAMETERS_NAME %s\
        -MAG_ZEROPOINT 25.0 -CATALOG_TYPE FITS_LDAC\
        " % (file_name, config_sex, file_name+'.ldac', param_sex)
        os.system(command)
        print('Executing command: %s\n' % command)
        
        
def run_scamp(ctext, catname):

    if os.path.exists('scampfilelist'):
        os.remove('scampfilelist')
        
                 
    scamp_list='scampfilelist'
    group_similar_files(scamp_list, common_text=ctext)
    cwd=os.getcwd()    
    config_scamp=cwd+'/config.scamp'
    
    
    scamp_command="scamp -c %s @%s" %(config_scamp, scamp_list)
    os.system(scamp_command)
    print('Executing command: %s\n' % scamp_command)       
        
   
   	



def run_swarp(ra, dec, pxscale, gain, ctext):
    
    if os.path.exists('wcsfiles'):
        os.remove('wcsfiles')

    wcs_files='wcsfiles'
    list_wcsfiles=group_similar_files(wcs_files, common_text=ctext)  
    cwd=os.getcwd()
    configfile=cwd+'/config.swarp'
    
    swarp_command="swarp -c %s @%s -SUBTRACT_BACK Y -RESAMPLE Y -RESAMPLE_DIR .\
    -COMBINE N -CENTER_TYPE MANUAL -CENTER %.3f, %.3f -PIXELSCALE_TYPE MANUAL\
    -PIXEL_SCALE %.3f -GAIN_DEFAULT %.3f" %(configfile,wcs_files,ra,dec,pxscale,gain)
        
        
    print("----------SwArP is running----------")
    print(swarp_command)
    os.system(swarp_command)
    
    
    
group_similar_files('list_object', '*.fits')  
run_astrometry('*.fits')  
list_wcsfiles=group_similar_files('list_wcsobject', 'wcs_*.fits')    
img_header=fits.open(list_wcsfiles[0])
data=img_header[0].data
header=img_header[0].header
w=WCS(header)
cd11=header['CD1_1']
cd21=header['CD2_1']
cd12=header['CD1_2']
cd22=header['CD2_2']
xscale=np.sqrt(cd11**2+cd21**2)
yscale=np.sqrt(cd12**2+cd22**2)
[ra,dec]=w.all_pix2world(data.shape[0]/2, data.shape[1]/2,1)
print('The Right Ascension of the Center of the field is:', ra)
print('The Declination of the Center of the field is:', dec)
pxscale=xscale*3600
print('The plate scale of the image is:', pxscale)
gain=1.6    
    

make_gaia_catalog(ra, dec, 1.0, 10, 22, catname='gaiacatalog')
run_sextractor(ctext='wcs_*.fits')
run_scamp(ctext='wcs_*.ldac', catname='gaiacatalog.ldac')
run_swarp(ra, dec, pxscale, gain, ctext='wcs_*.fits')

for text in ['*.xyls', '*.axy', '*.corr', '*.match', '*.new', '*.wcs', '*.solved', '*.rdls', '*.png', '*_resamp.weight.fits']:
    remove_similar_files(common_text=text)          

    
       
       
        
        
        
        
        

    
