# A python code for querying various catalogs
# Import the relevant packages and modules
# Author: Anirban Dutta
# Version: 1.0
# Date: 06/07/2019

# -------------------------------------------------------------------------------------------------------#

# Import the relevant packages and modules
import os
import glob
import shutil
from astropy.io import fits
from astroquery.vizier import Vizier
from astroquery.gaia import Gaia
import astropy.units as u
import astropy.coordinates as coord
from astropy.io import ascii
import warnings
from astropy.coordinates import SkyCoord
# -------------------------------------------------------------------------------------------------------#

# Ignoring warnings in output 

warnings.filterwarnings("ignore")

# -------------------------------------------------------------------------------------------------------#

"""
Functions to convert FITS files or astropy Tables to FITS_LDAC files and
vice versa.
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

'''
Functions to query GAIA and PANSTARSS around a given RA, DEC, search RADIUS, 
maximum magnitude and maximum sources

'''

def query_gaia(ra_deg, dec_deg, radius_deg, minmag=10, maxmag=20, maxsources=10000,catalogname='gaiacatalog'):
    job = Gaia.launch_job_async("SELECT * FROM gaiadr2.gaia_source AS g, gaiadr2.panstarrs1_best_neighbour \
    AS pbest, gaiadr2.panstarrs1_original_valid AS ps1 WHERE g.source_id = pbest.source_id AND \
    pbest.original_ext_source_id = ps1.obj_id AND CONTAINS(POINT('ICRS', g.ra, g.dec), \
        CIRCLE('ICRS', %.4f, %.4f, %.4f))=1 AND ps1.r_mean_psf_mag > %.2f AND ps1.r_mean_psf_mag \
        < %.2f AND pmra IS NOT NULL AND pmdec IS NOT NULL AND abs(pmdec) > 0 AND \
        abs(pmdec) < 40 AND abs(pmra)>0 AND abs(pmra) < 40 AND ps1.n_detections > 6 \
        AND pbest.number_of_mates=0 AND \ pbest.number_of_neighbours=1;"
        %(ra_deg, dec_deg, radius_deg, 10, 20))
	
	
    p = job.get_results()
    p['ra_errdeg'] = p['ra_error'] / 3.6e6
    p['dec_errdeg'] = p['dec_error'] / 3.6e6
    p['FLAGS'] = 0
    if os.path.exists(catalogname+'.ldac'):
        os.remove(catalogname+'.ldac')
    
    if os.path.exists('gaiacatalog.dat'):
        os.remove('gaiacatalog.dat')
    
    save_table_as_ldac(p, catalogname+'.ldac')

    

def query_PS1(ra_deg, dec_deg, radius_deg, maxmag=20, maxsources=5000,
catalogname='ps1catalog'):
    
    vquery=Vizier(columns=['objID','RAJ2000','e_RAJ2000','DEJ2000','e_DEJ2000','rmag','e_rmag'], 
    column_filters={"rmag":("<%f"%maxmag)}, row_limit=maxsources)

    field=coord.SkyCoord(ra=ra_deg, dec= dec_deg, unit=(u.deg,u.deg), frame='icrs')
    catalog=vquery.query_region(field,width=("%fd"%radius_deg),catalog="II/349/ps1")[0]
   
    if os.path.exists(catalogname+'.ldac'):
        os.remove(catalogname+'.ldac')
        
    save_table_as_ldac(catalog, catalogname+'.ldac')        
    
    if os.path.exists('ps1catalog.dat'):
        os.remove('ps1catalog.dat')
    ascii.write(catalog, 'ps1catalog.dat')
    
# -------------------------------------------------------------------------------------------------------#

'''
Function to remove files

'''

def remove_file(file_name):
    
    os.remove(file_name)
    
'''
Function to remove similar files
'''
    
def remove_similar_files(common_text):
    
    for file_name in glob.glob(common_text):
        remove_file(file_name)

'''
Function to make a list of the files

'''
def list_of_files(file_list, common_text):
        
    list_files=glob.glob(common_text)
    list_files.sort()
    with open(file_list, 'w') as f:
        for file_name in list_files:
            f.write(file_name+'\n')
            
            
    return list_files
    
'''
Function to convert a text list to a python list_files

'''    
def text_list_to_python_list(text_list):
    if os.path.exists(text_list):
        with open(text_list, 'r+') as f:
            python_list=f.read().split()
            
    return python_list
    
def python_list_to_text_list(python_list, text_list):
    with open(text_list, 'w') as f:
        for element in python_list:
            f.write(str(element)+'\n')
            
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
    list_files=list_of_files(file_list, common_text=ctext)

    
    remove_similar_files('*.cat')
    
    for file_name in list_files:
        command="sextractor %s -c %s -CATALOG_NAME %s -PARAMETERS_NAME %s\
        -MAG_ZEROPOINT 25.0 -CATALOG_TYPE FITS_LDAC\
        " % (file_name, config_sex, file_name+'.cat', param_sex)
        os.system(command)
        print('Executing command: %s\n' % command)
        
        

def run_scamp(ctext, ra, dec):

    if os.path.exists('scampfilelist'):
        os.remove('scampfilelist')
        
                 
    scamp_list='scampfilelist'
    list_of_files(scamp_list, common_text=ctext)
    cwd=os.getcwd()    
    config_swarp=cwd+'/config.scamp'
    
    
    scamp_command="scamp -c %s @%s" %(config_swarp, scamp_list)
    os.system(scamp_command)
    print('Executing command: %s\n' % scamp_command)
    
    
run_sextractor('*.fits')  
list_files=text_list_to_python_list('ImageList.list')
list_files.sort()

image=fits.open(list_files[0])
data=image[0].data
header=image[0].header
RA=image[0].header['TARRA']
DEC=image[0].header['TARDEC']

c= SkyCoord(ra=RA,dec=DEC,unit=(u.hourangle, u.deg))

ra= c.ra.deg
dec=c.dec.deg

 
#query_gaia(ra,dec,2)   
 
run_scamp('*.cat', ra, dec)
        
       
        
    
                        




















    
