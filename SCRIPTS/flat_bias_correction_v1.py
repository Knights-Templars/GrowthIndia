#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
    Created on Sun Feb 10 17:53:53 2019
    
    @author: hk
    """

from astropy.io import fits
import matplotlib.pyplot as plt
import glob
import os
from astropy.stats import sigma_clipped_stats
from astropy.stats import SigmaClip
import numpy as np
import numpy.ma as ma
import astropy.units as u
from astropy.io import ascii
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.stats import sigma_clipped_stats, sigma_clip
from astropy.time import Time
import astroscrappy
import argparse


parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("-d",help="Date of data to be reduces in YYYYMMDD format")
#parser.add_argument("-o",help="Object directory name",default='EMGW')
parser.add_argument("--box",help="Box size to search for catalog stars in arcminutes",type=float, default=30.0)
parser.add_argument("-m",help="Magnitude cut-offs of sources to be cross-matched against",type=float, default=19)
parser.add_argument("--anuradius",help="Magnitude cut-offs of sources to be cross-matched against",type=float, default=28/2.0)
parser.add_argument("--ra",help="RA of the target in decimal degrees",type=float)
parser.add_argument("--dec",help="DEC of the target in decimal degrees",type=float)
parser.add_argument("--a_min",help="Minimum value of aperture around the star",type=float, default=5)
parser.add_argument("--a_max",help="Minimum value of aperture around the star",type=float, default=15)
parser.add_argument("--process_files",help="Process already processed image",type=str ,default='N')
parser.add_argument("--cat",help="Process already processed image" ,default='ps1')

args = parser.parse_args()

date_dir=args.d
#obj_name=args.o
boxsize =args.box
maxmag=args.m
anuRadius =args.anuradius
ra=args.ra
dec=args.dec
a_min=args.a_min
a_max=args.a_max
process_files=args.process_files
cat=args.cat

#path to different directories

#path to different directories
default_path=("/home/anirban.dutta/Growth")
print(default_path)
basic_path=("/home/anirban.dutta/Growth")
flatFolder = os.path.join(basic_path+date_dir+"/flats/")
master_flat_folder=("/home/anirban.dutta/Growth/Mflats/")
biasFolder = (basic_path+date_dir+"/bias/")
print('Process Files flag', process_files)

'''
exclude_dir=['bias','panstarrs','observations.csv','bad_flats','GIT2019a', 'flat', 'reduced', '.DS_Store', 'dark', 'flats', 'focus']
object_list=[]
for direc in os.listdir(basic_path+date_dir):
    if direc not in exclude_dir:
        object_list.append(direc)

x_len=4108
y_len=4096


# Making masterbias
biasList = glob.glob(biasFolder+'*fits')
numBiasFiles = len(biasList)
print('Found %d bias files'%numBiasFiles)
biasImages = np.zeros((x_len, y_len, numBiasFiles))
for i in range(numBiasFiles):
    biasImages[:,:,i] = fits.open(biasList[i])[0].data
masterBias = np.median(biasImages, axis=2)   # gives problems with large images
print(masterBias)


# Making masterflats for different filters
print("flatfolder is ", flatFolder)
curpath = os.path.abspath('.')
if os.path.exists(flatFolder):
    print("Directory has flat files.")
    flag=False


else:
    flatFolder = master_flat_folder
    flag=True

os.chdir(flatFolder)
flatList = glob.glob('*.fits')



os.chdir(curpath)
flatlist_g=[]
flatlist_r=[]
flatlist_i=[]
flatlist_u=[]
flatlist_z=[]
flatlist_x=[]
flatlist=[flatlist_g,flatlist_r,flatlist_i,flatlist_u,flatlist_z,flatlist_x]
flatmap = {'u':0,'g':1, 'r':2, 'i':3, 'z':4, 'x':5}
for i, flatfilename in enumerate(flatList):
    filtername = fits.getval(flatFolder+flatfilename,'FILTER')
    try:
        flatlist[flatmap[filtername]].append(flatList[i])
    except KeyError:
        print("some weird unknown filter")
x=[]
for j, list in enumerate(flatlist):
    if not list:
        x.append(j)
print(x)
if x:
    os.chdir(master_flat_folder)
    flatList_1 = glob.glob('*.fits')
    os.chdir(curpath)
    
    for i, flatfile in enumerate(flatList_1):
        filtername = fits.getval(master_flat_folder+flatfile,'FILTER')
        if (flatmap[str(filtername)]) in x:
            try:
                flatlist[flatmap[filtername]].append(flatList_1[i])
            except KeyError:
                print("some weird unknown filter")
print(flatlist)

# Function to combine flat frames
def flatcombine(image_list, masterBias, flatFolder):
    """
        Take a list of flats and combine them
        """
    num_files = len(image_list)
    flatImages = np.empty((num_files, x_len, y_len))
    
    for c_files, this_file in enumerate(image_list):
        data = fits.getdata(flatFolder+this_file) - masterBias
        data_med = np.median(data)
        flatImages[c_files,:,:] = data / data_med
    return np.median(flatImages, axis=0)

masterFlat = np.zeros([len(flatmap), x_len, y_len])
# create a flat for each filter
for count, filterflatlist in enumerate(flatlist):
    if count in x:
        flatFolder = master_flat_folder
    else:
        if not flag:
            flatFolder=  os.path.join(basic_path+date_dir+"/flats/")
    if (len(filterflatlist)>0):
        masterFlat[count, :, :] = flatcombine(filterflatlist, masterBias, flatFolder)
    print(np.median(masterFlat[count]))


for obj_name in object_list:
    sciFolder =  os.path.join(basic_path,date_dir,obj_name+'/')
    procFolder = os.path.join(basic_path,date_dir,obj_name+"/reduced/")
    print(sciFolder)


    print("Path of science directory is : ",sciFolder)
    curpath = os.path.abspath('.')
    os.chdir(sciFolder)
    sciList = glob.glob('*wcs.fits')
    os.chdir(curpath)
    numSciFiles = len(sciList)
    if not os.path.exists(sciFolder) :
        print("Directory is empty. Please check the path you have specified.")
    else:
        print('Found %d science files'%numSciFiles)
        print(sciList)
    scilist_g=[]
    scilist_r=[]
    scilist_i=[]
    scilist_u=[]
    scilist_z=[]
    scilist_x=[]
    scilist=[scilist_g,scilist_r,scilist_i,scilist_u,scilist_z,scilist_x]
    scimap = {'u':0,'g':1, 'r':2, 'i':3, 'z':4,'x':5}
    for i, scifilename in enumerate(sciList):
        filter_name = fits.getval(sciFolder+scifilename,'FILTER')
        print(filter_name)
        try:
            scilist[scimap[filter_name]].append(sciList[i])
        except KeyError:
            print("some weird unknown filter")
    
    detectorGain = 1.6 #in e-/ADU
    readNoise = 14.0 #in electrons
    saturLevel = 150000 #in electrons
    if os.path.exists(procFolder):
        print("reduced folder already exists.")
    else:
        os.system("cp -r "+default_path+"reduced"+" "+basic_path+date_dir+"/"+obj_name+"/")
        procFolder = os.path.join(basic_path,date_dir,obj_name+"/reduced/")
    print(procFolder)
    
    import subprocess
    #from subprocess import call
    curpath = os.path.abspath('/mnt/growth/')
    os.chdir(curpath)
    counter=0
    #Dan Perley's astrometry python code
    autoastrometry_script = curpath+'/autoastrometry.py'
    os.chdir(procFolder)
    red_list= glob.glob('*.proc.cr.fits')
    for k in range(len(scilist)):
        if ((np.median(masterFlat[k]))> 0):
            print("Entered loop 2")
            for i in range(len(scilist[k])):
                print(scilist[k][i])
                img=scilist[k][i] +'.proc.cr.fits'
                print('image name =', img)
				
                if img in red_list and process_files == 'N':
                    print('whoooohoooo Already readuced. use process_files keyword to reduce it again')
                   
                    processed=True
                    continue
                else:
                    
                    rawHDU = fits.open(sciFolder+scilist[k][i])[0]
                    rawData = rawHDU.data
                    rawHeader = rawHDU.header
                    
                    #Correct for the bias and flats here
                    procData = (rawData - masterBias)/masterFlat[k]
                
                    #           print(rawData[1964][2169],'///',masterBias[1965][2170],'///',masterFlat[k,2170,1965],'///',procData[1964][2169])
                    procHDU = fits.PrimaryHDU(procData)
                    #Write the reduced frame to disk, propagating the original header of the raw image
                    # procHDU.writeto(procFolder+scilist[k][i]+'.proc.fits', overwrite=True)
                    crmask, cleanArray = astroscrappy.detect_cosmics(procData,gain=detectorGain,readnoise=readNoise, satlevel=saturLevel)
                    
                    print('Number of affected pixels is %d for file %s'%(np.sum(crmask), scilist[k][i])) #print number of affected pixels
                     
                    #The returned clean array is in units of electrons -- have to divide it by the gain to recover ADUs
                    procData_cr = cleanArray / detectorGain
                    crCleanHDU = fits.PrimaryHDU(procData_cr)
                    crCleanHDU.header = rawHeader
		    crCleanHDU.header.remove('BZERO')
		    crCleanHDU.header.remove('BSCALE')
                    crCleanHDU.writeto(procFolder+ scilist[k][i] +'.proc.cr.fits', overwrite=True)
                   
                    
                    import subprocess
                    #for k in range(len(scilist)):
                    #if ((np.median(masterFlat[k]))> 0):
                    # for i in range(len(scilist[k])
                
                #fig, ax = plt.subplots(1, figsize=(8,8))
                    procData = fits.open(scilist[k][i]+'.proc.cr.fits')[0].data
                    mean, median, std = sigma_clipped_stats(procData)
                #plt.imshow(procData, vmin = median - 3*std, vmax = median + 3*std)
                
                #overlay the asterisms on the plot
                #for p in patch_list:
                #    ax.add_patch(p)
                #for t in artist_list:
                #    ax.add_artist(t)
                
                #plt.colorbar()
                #plt.show()
                
                    print(k, i, scilist[k][i])
                    imageName= scilist[k][i]+'.proc.cr.fits'
                
                    print (imageName)
                    os.system("pwd")
                    data, header = fits.getdata(imageName, header=True)
                    filtername=header['FILTER']

                    if args.ra and args.dec is not None:
                        ra = args.ra
                        dec = args.dec
                    else:
                        coords=SkyCoord(ra=header['TARRA'], dec=header['TARDEC'], unit=(u.hour,u.deg))
                        ra=coords.ra.degree
                        dec=coords.dec.degree
                    date_t=(header["DATE-OBS"]).split('T')
                    date= str(date_t[0])+ ' ' + str(date_t[1])
                    print(date)
                #Compute some image statistics for scaling the image plot
                    mean, median, sigma = sigma_clipped_stats(data)
                #plot the image with some reasonable scale
                #plt.figure(figsize=(10,10))
                #plt.imshow(data, vmin=median-3*sigma, vmax=median+3*sigma)
                #plt.show()
                    w = WCS(header)
                #Get the RA and Dec of the center of the image
                    [raImage, decImage] = w.all_pix2world(data.shape[0]/2, data.shape[1]/2, 1)
                    print([raImage, decImage])
                #Set the box size to search for catalog stars
                # arcminutes
                #Magnitude cut-offs of sources to be cross-matched against
                    from astroquery.vizier import Vizier
                #Vizier.VIZIER_SERVER = 'vizier.ast.cam.ac.uk'
                    if cat!='ps1':
                        catNum = 'V/147'
                    else:
                        catNum = 'II/349'#This is the catalog number of PS1 in Vizier
                    print('\nQuerying Vizier %s around RA %.4f, Dec %.4f with a radius of %.4f arcmin'%(catNum, raImage, decImage, boxsize))
                #
                #You can set the filters for the individual columns (magnitude range, number of detections) inside the Vizier query
                    v = Vizier(columns=['*'], column_filters={filtername+"mag":"<%.2f"%maxmag, "Nd":">6", "e_"+filtername+"mag":"<<1.086/3"}, row_limit=-1)
                    Q = v.query_region(SkyCoord(ra = raImage, dec = decImage, unit = (u.deg, u.deg)), radius = str(boxsize)+'m', catalog=catNum, cache=False)
                #query vizier around (ra, dec) with a radius of boxsize
                    print(Q[0])
                #except:
                #    print('I cannnot reach the Vizier database. Is the internet working?')
                    if cat!='ps1':
                        ps1_imCoords = w.all_world2pix(Q[0]['RA_ICRS'], Q[0]['DE_ICRS'], 1)
                    else:
                        ps1_imCoords = w.all_world2pix(Q[0]['RAJ2000'], Q[0]['DEJ2000'], 1)
                #Another round of filtering where we reject sources close to the edges
                    good_cat_stars = Q[0][np.where((ps1_imCoords[0] > 500) & (ps1_imCoords[0] < 3500) & (ps1_imCoords[1] > 500) & (ps1_imCoords[1] < 3500))]
                    print(good_cat_stars)
                    if cat!='ps1':
                        ps1_imCoords = w.all_world2pix(good_cat_stars['RA_ICRS'],good_cat_stars['DE_ICRS'], 1)
                    else:
                        ps1_imCoords = w.all_world2pix(good_cat_stars['RAJ2000'],good_cat_stars['DEJ2000'], 1)

                    configFile = 'photomCat.sex'
                    catalogName = imageName+'.cat'
                    paramName = 'photomCat.param'
                    try:
                        rval = ('sextractor -c '+ configFile+' '+ imageName+' '+ '-CATALOG_NAME ' +catalogName+  ' -PARAMETERS_NAME '+ paramName)
		        rval=os.system(rval)
                    except:
                        print('Could not run sextractor with exit error %s'%err)
                    import astromatic_wrapper as aw
                #This is a python wrapper for reading LDAC files produced by SExtractor
                    sourceTable = aw.utils.ldac.get_table_from_ldac(catalogName)
                #Let's look at the contents of the table
                    print(sourceTable.colnames)
                    print(sourceTable)
                #filter on the sources to select the ones satisfying our criteria
                    cleanSources = sourceTable[(sourceTable['FLAGS']==0) & (sourceTable['FWHM_WORLD'] < 2) & (sourceTable['XWIN_IMAGE']<3500) & (sourceTable['XWIN_IMAGE']>500) &(sourceTable['YWIN_IMAGE']<3500) &(sourceTable['YWIN_IMAGE']>500)]
                    print(cleanSources.colnames)
                    print(cleanSources)
                #fig = plt.figure(figsize=(10,10))
                #ax = fig.gca()
                #plt.imshow(data, vmin=median-3*sigma, vmax=median+3*sigma)
                #plotting circles on top of all detected sources
                #circles = [plt.Circle((source['XWIN_IMAGE'], source['YWIN_IMAGE']), radius = 5, edgecolor='r', facecolor='None') for source in cleanSources]
                #for c in circles:
                #    ax.add_artist(c)
                #plt.show()
                    sourceCatCoords = SkyCoord(ra=cleanSources['ALPHAWIN_J2000'], dec=cleanSources['DELTAWIN_J2000'], frame='icrs', unit='degree')
                    if cat!='ps1':
                        ps1CatCoords = SkyCoord(ra=good_cat_stars['RA_ICRS'], dec=good_cat_stars['DE_ICRS'], frame='icrs', unit='degree')
                    else:
                        ps1CatCoords = SkyCoord(ra=good_cat_stars['RAJ2000'], dec=good_cat_stars['DEJ2000'], frame='icrs', unit='degree')
                   
                  #Now cross match sources
                #Set the cross-match distance threshold to 0.6 arcsec, or just about one pixel
                    photoDistThresh = 1.0
                    idx_image, idx_ps1, d2d, d3d = ps1CatCoords.search_around_sky(sourceCatCoords, photoDistThresh*u.arcsec)
                #idx_image are indexes into sourceCatCoords for the matched sources, while idx_ps1 are indexes into ps1CatCoords for the matched sources
                    print('Found %d good cross-matches'%len(idx_image))
                #file = open("test.txt","w")
                #print(idx_imag)
                #for i in (idx_ps1):
                # print(str(ps1CatCoords[i]))
                #file.close()
                    aperture_diameter = np.arange(a_min,a_max)
                #For each aperture, we are going to compute the magniutde difference between the largest pixel aperture and that specific aperture for every source in cross-matched catalog
                    try:

    		        magDiff = np.ma.zeros((len(aperture_diameter), len(idx_image)))
                        sig_clip_1=SigmaClip(sigma=1)
                        sig_clip_2=SigmaClip(sigma=2) 
                        for j in range(len(aperture_diameter)):
                            magDiff[j] = sigma_clip(cleanSources['MAG_APER'][:,9][idx_image] - cleanSources['MAG_APER'][:,j][idx_image])
                            zeroPoints = []
                            for i in range(len(aperture_diameter)):
                        #Array of differences between the catalog and instrumental magnitudes
                                offsets = ma.array(good_cat_stars[filtername+'mag'][idx_ps1] - cleanSources['MAG_APER'][:,i][idx_image])
                                moffsets_1= sig_clip_1(offsets)
                                moffsets_2= sig_clip_2(offsets)
                        
                        #plt.plot(offsets)
                        
                        #Compute sigma clipped statistics
                                if(len(moffsets_1[moffsets_1.mask==False]))>=15:
                                    zero_mean, zero_med, zero_std = sigma_clipped_stats(offsets,sigma=1)
                                else:
                                    if (len(moffsets_2[moffsets_2.mask==False]))>=15:
                                        zero_mean, zero_med, zero_std = sigma_clipped_stats(offsets,sigma=2)
                                    else:
                                        zero_mean, zero_med, zero_std = sigma_clipped_stats(offsets)
                                zeroDict = {'diameter': aperture_diameter[i], 'zp_mean': zero_mean, 'zp_median': zero_med, 'zp_std': zero_std}
                                zeroPoints.append(zeroDict)
                        print((zeroDict))
                        medFWHM = np.median(cleanSources['FWHM_WORLD'][idx_image])
                        print('The median FWHM in the image is %.2f arcsec'%(medFWHM*3600))
                    
                        MJD=np.round(Time(date).mjd,3)
                        print(date,MJD)
                        header['ZP']=(zeroDict['zp_median'], 'Zero point obtained after reduction')
                        header['ZP_err']=(zeroDict['zp_std'], 'Error in zero point obtained after reduction')
                        header['N_Source']=( len(idx_image), 'No. of sources used to find zero point')
                        header['med_FWHM']=( medFWHM*3600, 'median FWHM in the image/seeing')
                        fits.writeto(imageName, data, header, overwrite=True)
                
                        f=open('reduction_table.txt','a+')
                        f.write(' \n ')
                        if (counter==0):
                            f.write('Object\tMJD/t_0\tFilter\n')
                        f.write('{}\t{}\t{}\t{}\n'.format(obj_name,MJD,filtername, medFWHM*3600))
                        f.close()
		    except:
			continue

'''
