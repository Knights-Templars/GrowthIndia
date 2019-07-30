# Imprrt required libraries

import os
import re
import sys
import glob
import shutil
from pyraf import iraf
from astropy.io import fits


# Global Variables

filters_headname=['7BesU', '6BesB', '5BesV', '4BesR', '3BesI']
filters=[x[-1] for x in filters_headname]

OBJECT_NAME='2017hpa'

# Telescope CCD specifications

read_noise=4.87
gain=1.22
data_max=55000

#image header keyword

Date='DATE-OBS'
Filter='FILTER'
OBJECT='OBJECT'
DIR_SN='/home/anirban.dutta/NEW/'
DIR_PHOT=DIR_SN+"/SN_Photometry/"
DIR_PG=DIR_SN+"PG_Photometry/"

# Functions for file handling

# function for removing files: argument: file_name

def remove_file(file_name):
	
	try:
		os.remove(file_name)
	except OSError:
	    pass

# Function for removing files having similar names: argument: common_text

def remove_similar_files(common_text):
	
	for residual_file in glob.glob(common_text):
		remove_file(residual_file)
		
		
for text in ['*b_flat*', 'list_*']:
    remove_similar_files(common_text=text)

# Function to group similar files: arguments: text_list, common_text, exceptions

def group_similar_files(text_list, common_text, exceptions=''):
	list_files=glob.glob(common_text)
	if exceptions !='':
		list_exceptions=exceptions.split(',')
		for file_name in glob.glob(common_text):
			for text in list_exceptions:
				test=re.search(text, file_name)
				if test:
					try:
						list_files.remove(file_name)
					except ValueError:
						pass
	list_files.sort()
	if len(text_list) !=0:
		with open(text_list, 'w') as f:
			for file_name	in list_files:
				f.write(file_name + '\n')

	return list_files


def text_list_to_python_list(text_list):
	if os.path.exists(text_list):
		with open(text_list, 'r+') as f:
			python_list=f.read().split()

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
                
            
            

# Load IRAF Packages:

iraf.noao(_doprint=0)
iraf.imred(_doprint=0)
iraf.ccdred(_doprint=0)
iraf.crutil(_doprint=0)
iraf.images(_doprint=0)
iraf.ccdred.instrument='ccddb$kpno/direct.dat'

def hedit(textlist_files, field_name, value, add_keyword ='no'):

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
    
    
def zero_combine(textlist_bias, master_bias):
    
    task=iraf.noao.imred.ccdred.zerocombine
    task.unlearn()
    
    task.combine='median'
    task.reject='ccdclip'
    task.ccdtype=''
    task.process='no'
    task.delete='no'
    task.clobber='no'
    task.scale='none'
    task.statsec=''
    task.nlow=0
    task.nhigh=1
    task.nkeep=1
    task.mclip='yes'
    task.lsigma=3.
    task.hsigma=3.
    task.rdnoise=float(read_noise)
    task.gain=float(gain)
    task.snoise=0
    task.pclip=-0.5
    task.blank=0
    task.mode='ql'
    
    remove_file(master_bias)
    task(input='@'+textlist_bias, output=master_bias)

def bias_subtract(textlist_tbs, master_bias='mbias.fits', prefix_str='b_'):

    task=iraf.noao.imred.ccdred.ccdproc
    task.unlearn()
    
    task.ccdtype=''
    task.max_cac=0
    task.noproc='no'
    task.fixpix='no'
    task.oversca='no'
    task.trim='no'
    task.zerocor='yes'
    task.darkcor='no'
    task.flatcor='no'
    task.illumco='no'
    task.fringec='no'
    task.readcor='no'
    task.scancor='no'
    task.readaxi='line'
    task.fixfile=''
    task.biassec=''
    task.trimsec=''
    task.zero=master_bias
    task.dark=''
    task.flat=''
    task.fringe=''
    task.minrepl=1
    task.scantyp='shortscan'
    task.nscan=1
    task.interac='no'
    task.functio='legendre'
    task.order=1
    task.sample='*'
    task.naverag=1
    task.niterat=1
    task.low_rej=3
    task.high_re=3
    task.grow=0
    task.mode='ql'
    
    output_filename=prefix_str+'//'+'@'+textlist_tbs
    remove_file(output_filename)
    
    task(images='@'+textlist_tbs, output=output_filename)
    
def image_statistics(ctext, text_list='list_bsflat'):

    textlist_bsflat='list_bsflat'
    group_similar_files(textlist_bsflat, common_text=ctext)
    
    task=iraf.images.imutil.imstatistics
    task.unlearn()
    
    task.lower='INDEF'
    task.upper='INDEF'
    task.nclip=0
    task.lsigma=3
    task.usigma=3
    task.binwidt=0.1
    task.cache='no'
    task.mode='ql'
    
    task(format='no', images='@'+textlist_bsflat, field='mean', Stdout=text_list+'_mean')
    task(format='yes', images='@'+textlist_bsflat, field='image,mean,midpt,stddev,min,max', Stdout=text_list+'_stat')
    
    list_mean=text_list_to_python_list(text_list+'_mean')
    
    return list_mean
    
def flat_normalize(ctext, textlist_mean='list_bsflat', prefix_str='n'):

    list_bsflat=group_similar_files('', common_text=ctext)
    list_bsflat_mean=image_statistics(ctext, textlist_mean)
    
    task=iraf.images.imutil.imarith
    task.unlearn()
    
    for index in range(0, len(list_bsflat)):
        output_filename=prefix_str+list_bsflat[index]
        remove_file(output_filename)
        task(operand1=list_bsflat[index], op='/', operand2=float(list_bsflat_mean[index]), result=output_filename)
        

def flat_combine(ctext):
    
    band=ctext[-7]
    textlist_nbflat='list_nbflat'
    group_similar_files(textlist_nbflat, common_text=ctext)
    
    task=iraf.noao.imred.ccdred.flatcombine
    task.unlearn()
    
    task.combine='median'
    task.reject='ccdclip'
    task.ccdtype=''
    task.process='no'
    task.subsets='no'
    task.delete='no'
    task.clobber='no'
    task.scale='mode'
    task.statsec=''
    task.nlow=1
    task.nhigh=1
    task.nkeep=1
    task.mclip='yes'
    task.lsigma=3
    task.hsigma=3
    task.rdnoise=float(read_noise)
    task.gain=float(gain)
    task.snoise=0
    task.pclip=-0.5
    task.blank=1
    task.mode='ql'
    
    output_filename='mflat'+band+'.fits'
    remove_file(output_filename)
    task(input='@'+textlist_nbflat, output=output_filename)

def flat_correction(ctext, master_flat, exception='flat', prefix_str='f'):
    
    list_bsobject=group_similar_files('',common_text=ctext, exceptions=exception)
    
    task=iraf.images.imutil.imarith
    task.unlearn()
    
    for image in list_bsobject:
        output_filename=prefix_str+image
        remove_file(output_filename)
        task(operand1=image, op='/', operand2=master_flat, result=output_filename)
        
        
def crmedian(ctext, prefix_str='c'):
    
    list_cosmic=group_similar_files('', common_text=ctext)
    
    task=iraf.noao.imred.crutil.crmedian
    task.unlearn()
    
    task.crmask=''
    task.median=''
    task.sigma=''
    task.residua=''
    task.var0=0
    task.var1=0
    task.var2=0
    task.lsigma=10
    task.hsigma=3
    task.ncmed=5
    task.nlmed=5
    task.ncsig=25
    task.nlsig=25
    task.mode='ql'
    
    for image in list_cosmic:
        output_filename=prefix_str+image
        remove_file(output_filename)
        task(input=image, output=output_filename)
    
    
remove_prev_files=True
if remove_prev_files:
    for text in ['tmp', '*b_*', 'mbias', 'mflat*', 'list_']:
        remove_similar_files(common_text=text)
        
# Group Similar types of .fits files

group_similar_files('list_bias', 'Bias*.fits')
group_similar_files('list_flat', 'flat*.fits')
group_similar_files('list_object', OBJECT+'-*.fits')
group_similar_files('list_standard', 'PG*.fits')

for band in filters:
    group_similar_files('list_'+band.upper(), '*-'+band.upper()+'*.fits')     
    
list_list_phot=['list_bias', 'list_flat', 'list_object', 'list_standard']
list_lists_to_list(list_list_phot[1:], 'list_phot')
           
    
hedit('list_bias', OBJECT, 'FLAT')
hedit('list_flat', OBJECT, 'BIAS')
hedit('list_object', OBJECT, str(OBJECT_NAME))

for index, band in enumerate(filters):
    hedit('list_'+band.upper(), Filter, filters_headname[index])
    

# combine the bias to make a master bias and apply bias correction

zero_combine(textlist_bias='list_bias', master_bias='mbias.fits')
bias_subtract(textlist_tbs='list_phot')

flat_normalize(ctext='b_flat*.fits')

for band in filters:
    flat_combine('nb_flat_'+band.upper()+'*.fits')
    
for band in filters:
    flat_correction(ctext='b_*_'+band.upper()+'*.fits', master_flat='mflat'+band.upper()+'.fits')
    
crmedian(ctext='fb_*.fits')

if os.path.exists(DIR_PHOT):
    os.remove(DIR_PHOT)

os.mkdir(DIR_PHOT)
    


for file_name in group_similar_files('', common_text='cfb_*.fits'):
    header=fits.getheader(file_name, ext=0)
    date_obs=header[Date]
    if os.path.isfile(DIR_PHOT+date_obs+'-'+file_name):
        os.remove(DIR_PHOT+date_obs+'-'+file_name)
    shutil.copy(file_name, DIR_PHOT+date_obs+'-'+file_name)
    

    
    
       


