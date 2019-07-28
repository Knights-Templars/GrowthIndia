# A Python Code for performing Subtraction of Ps1 images from Growth-India Telescope images having WCS and aligned by SwArP 
# Author - Anirban Dutta, Avinash Singh
# Version - 1.0
# Date: 16/07/2019

#----------------------------------------------------------------------#
# Import relevant packages and modules

import os
import glob
import shutil
import re
import numpy as np
from astropy.io import ascii
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.wcs import WCS
import warnings


from pyraf import iraf

#----------------------------------------------------------------------#

# ---------------------------------------------------------------------#

# Ignoring warnings in output 

warnings.filterwarnings("ignore")

# ---------------------------------------------------------------------#


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
        
# Function for grouping similar files
       
def group_similar_files(text_list, common_text, exceptions=''):
    
	list_files=glob.glob(common_text)
	if exceptions !='':
		list_files=filter(lambda z: not re.search(exceptions, z), list_files)

	list_files.sort()
	if len(text_list) !=0:
		with open(text_list, 'w') as f:
			for file_name in list_files:
				f.write(file_name+'\n')
                
	return list_files    
    
# Function to convert a text list to a python list

def text_list_to_python_list(text_list):
    if os.path.exists(text_list):
        with open(text_list, 'r+') as f:
            python_list=f.read.split()
            
    return python_list
    
# Function to convert a python list to text list

def python_list_to_text_list(python_list, text_list):
    with open(text_list, 'w') as f:
        for element in python_list:
            f.write(str(element)+'\n')
            
# Function to group file names from a list to a single text file           
            
def list_lists_to_list(list_lists, text_list):
    list_name=[]
    for file_name in list_lists:
        with open(file_name, 'r') as f:
            file_list=f.read().split()
            for element in file_list:
                list_name.append(element)
    python_list_to_text_list(list_name, text_list)
    
    return list_name    
    
#----------------------------------------------------------------------#

# Get the current working directory


#----------------------------------------------------------------------#

# Remove similar files from previous run

for text in ['*.cat', '*.coo', 'conv_*.fits', 'sub_*.fits', 'ker*.fits', 'ref_*.fits', '*.list']:
	remove_similar_files(common_text=text)

#----------------------------------------------------------------------#

# Function to run sextractor on the first file in the science images list and make a coordinate file to be used by the task PsfMatch

def run_sextractor_coord(sci_image_name):
	
	cwd= os.getcwd()
	config_sex=cwd+'/config.sex'
	param_sex=cwd+'/sextractor.param'

	print("Running SeXtRaCtOr")
	if os.path.exists('*.cat'):
		os.remove('*.cat')

	catname=sci_image_name+'.cat'
	command="sextractor %s -c %s -CATALOG_NAME %s -MAG_ZEROPOINT 25.0 -SATUR_LEVEL 55000" % (sci_image_name, config_sex, catname)
	os.system(command)
	print('Executing command: %s\n' % command)
	print('Success!')
	
	sci_cat=ascii.read(catname)
	clean_sci_sources=sci_cat[(sci_cat['FLAGS']==0) & (sci_cat['SNR_WIN'] > 100) & (sci_cat['SNR_WIN'] < 500)]
	coord_sci=clean_sci_sources['XWIN_IMAGE', 'YWIN_IMAGE']
	ascii.write(coord_sci, 'stars.coo')
	fwhm_sci=np.median(clean_sci_sources['FWHM_IMAGE'])
	elongation_sci=np.median(clean_sci_sources['ELONGATION'])
	ratio=1.0/elongation_sci
	position_angle=np.median(clean_sci_sources['THETAWIN_IMAGE'])

	if position_angle<0:
		corr_PA=180+position_angle
	else:
		corr_PA=position_angle
	

#----------------------------------------------------------------------#

# Function to use the IRAF task Psfmatch

def psf_match(input_image, reference_image, reference_stars, prefix_str='ref'):
	
	task=iraf.immatch.psfmatch
	task.unlearn()

	task.psfdata=reference_stars
	task.kernel=''
	task.convolu='image'
	task.dnx=31
	task.dny=31
	task.pnx=15
	task.pny=15
	task.center='yes'
	task.backgro='median'
	task.lorejec='INDEF'
	task.hirejec='INDEF'
	task.apodize=0.0
	task.fluxratio='INDEF'
	task.filter='cosbell'
	task.sx1='INDEF'
	task.sx2='INDEF'
	task.sy1='INDEF'
	task.sy2='INDEF'
	task.radsym='yes'
	task.thresho=0.2
	task.normfactor=1.0
	task.boundary='nearest'
	task.constant=0.0
	task.interac='no'
	task.verbose='yes'
	task.graphic='stdgraph'
	task.display='stdimage'
	task.gcomman=''
	task.icomman=''
	task.mode='ql'
	
	output_filename='conv_'+reference_image
	remove_file(output_filename)
	
	task(input=input_image, referenc=reference_image, psfdata=reference_stars, kernel='', output=output_filename)
	
	return output_filename

#------------------------------------


def run_sextractor_multifile(reference_image, ctext_sci, ctext_ref):
	
	cwd= os.getcwd()
	config_sex=cwd+'/config.sex'
	
	file_list='Imagelist.list'
	list_files=group_similar_files(file_list, common_text=ctext_sci)
	sci_image_name=list_files[0]
	run_sextractor_coord(sci_image_name)
	
	
	for file_name in list_files:
		
		catname=file_name+'.cat'
		cmd_sex_sci="sextractor -c %s %s -CATALOG_NAME %s -MAG_ZEROPOINT 25.0 -CATALOG_TYPE=ASCII_HEAD" % (config_sex, file_name, catname) 
		os.system(cmd_sex_sci)
		print('Executing command: %s\n' % cmd_sex_sci)
		print('Success!')
		psf_match(reference_image, file_name, reference_stars='stars.coo')
		for text in ['ker*.fits']:
			remove_similar_files(common_text=text)
		
		
	file_list='Reference_convolved.list'
	list_files=group_similar_files(file_list, common_text=ctext_ref)
	
	for file_name in list_files:
		
		catname_ref=file_name+'.cat'
		cmd_sex_ref="sextractor -c %s %s -CATALOG_NAME %s -MAG_ZEROPOINT 25.0 -CATALOG_TYPE=ASCII_HEAD" % (config_sex, file_name, catname_ref)
		os.system(cmd_sex_ref)
		print('Executing command: %s\n' % cmd_sex_ref)
		
		
run_sextractor_multifile('reference_*.fits', ctext_sci='wcs_*.fits', ctext_ref='conv_*.fits')
		

def flux_calibration(ctext_ref, ctext_sci):
	
	file_list_sci='ScienceCatalogs.list'
	list_files_sci=group_similar_files(file_list_sci, common_text=ctext_sci)
	
	file_list_ref='ConvolvedReferenceCatalogs.list'
	list_files_ref=group_similar_files(file_list_ref, common_text=ctext_ref)
	
	for file_name in list_files_sci:
		cat_sci=ascii.read(file_name)
		print("We are reading the science catalog:", file_name)
		cat_ref=ascii.read('conv_'+file_name)
		print("We are reading the reference catalog:", 'conv_'+file_name)
		c1 = SkyCoord(ra=cat_ref['X_WORLD'], dec=cat_ref['Y_WORLD'],unit='degree')
		c2 = SkyCoord(ra=cat_sci['X_WORLD'], dec=cat_sci['Y_WORLD'],unit='degree')
		idx_ref, idx_sci, d2d, d3d = c2.search_around_sky(c1, 1.0*u.arcsec)
		print('Number of cross-matched sources is %d'%(len(d2d)))
		flux_sci = cat_sci['FLUX_AUTO']
		flux_ref = cat_ref['FLUX_AUTO']
		scale = np.median(flux_sci[idx_sci]/ flux_ref[idx_ref])
		print("The scale is:", scale)
		
		sci_image=file_name[:-4]
		print("We are reading the convolved science image:", sci_image)
		sci_image_name=fits.open(sci_image)
		data_sci=sci_image_name[0].data
		hdr_sci=sci_image_name[0].header
		
		ref_image='conv_'+file_name[:-4]
		print("We are reading the convolved reference image:", ref_image)
		ref_image_name=fits.open(ref_image)
		ref_conv=ref_image_name[0].data
		
		image_sub=data_sci-scale*ref_conv
		hdu_image_sub=fits.PrimaryHDU(image_sub, hdr_sci)
		hdu_image_sub.writeto('sub_'+file_name[:-4])
		
		
flux_calibration(ctext_ref='conv_wcs_*.cat', ctext_sci='wcs_*.cat')
