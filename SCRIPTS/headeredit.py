import os
import numpy
import shutil
import glob
from astropy.io import fits
from pyraf import iraf
import subprocess
import numpy
from astropy.stats import sigma_clipped_stats


# Telescope CCD Specifications:
Telescope_name= 'HCT'
CCD_name='HFOSC_old'
read_noise=4.87
gain=1.22
data_max=55000
OBJECT='OBJECT'
OBJECT_NAME='PG0918+029'
Right_Ascension='RA'
Declination='DEC'
RA='09:21:28'
DEC='02:46:03'

iraf.images(_doprint=0)

cwd=os.getcwd()
config_dir=cwd+"/CONFIG/"
DIR_aligned=cwd+"/SN_ALIGNED/"

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
                    
# function to edit the header

def hedit(textlist_files, field_name, value, add_keyword='yes'):

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
    
for text in ['list_']:
    remove_similar_files(common_text=text)       
            
group_similar_files('list_object', 'PG*.fits')

hedit('list_object', Right_Ascension, str(RA), add_keyword='yes')
hedit('list_object', Declination, str(DEC), add_keyword='yes')
hedit('list_object', OBJECT, str(OBJECT_NAME), add_keyword='yes')
    
