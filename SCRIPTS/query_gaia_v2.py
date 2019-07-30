import os
from astropy.io import fits
from astropy.io import ascii
import warnings
import glob
from astropy.coordinates import SkyCoord
import astropy.units as u
warnings.filterwarnings("ignore")
from astroquery.gaia import Gaia
from astropy.wcs import WCS


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


def make_gaia_catalog(ra, dec, catalog_box_size, catalog_min_mag, catalog_max_mag, catname, writetext = True, writeldac = True):
	job = Gaia.launch_job_async("SELECT * FROM gaiadr2.gaia_source AS g, gaiadr2.panstarrs1_best_neighbour AS pbest, gaiadr2.panstarrs1_original_valid AS ps1 WHERE g.source_id = pbest.source_id AND pbest.original_ext_source_id = ps1.obj_id AND CONTAINS(POINT('ICRS', g.ra, g.dec), CIRCLE('ICRS', %.4f, %.4f, %.4f))=1 AND ps1.r_mean_psf_mag > %.2f AND ps1.r_mean_psf_mag < %.2f AND pmra IS NOT NULL AND pmdec IS NOT NULL AND abs(pmdec) > 0 AND abs(pmdec) < 40 AND abs(pmra)>0 AND abs(pmra) < 40 AND ps1.n_detections > 6 AND pbest.number_of_mates=0 AND pbest.number_of_neighbours=1;"%(ra, dec, catalog_box_size, catalog_min_mag, catalog_max_mag), dump_to_file = False)
	
	p = job.get_results()
	
	# convert RA and DEC errors from mas(milli arc second) to degrees

	p['ra_errdeg'] = p['ra_error'] / 3.6e6
	p['dec_errdeg'] = p['dec_error'] / 3.6e6
	p['FLAGS'] = 0
	
	p.remove_columns(['astrometric_n_obs_al', 'astrometric_n_obs_ac', 'astrometric_n_good_obs_al', 'astrometric_n_bad_obs_al', 'astrometric_gof_al', 'astrometric_chi2_al', 'astrometric_excess_noise', 'astrometric_excess_noise_sig', 'astrometric_params_solved', 'astrometric_primary_flag', 'astrometric_weight_al', 'astrometric_pseudo_colour', 'astrometric_pseudo_colour_error', 'mean_varpi_factor_al', 'astrometric_matched_observations', 'visibility_periods_used', 'astrometric_sigma5d_max', 'frame_rotator_object_type', 'matched_observations', 'duplicated_source', 'phot_g_n_obs', 'phot_g_mean_flux', 'phot_g_mean_flux_error', 'phot_g_mean_flux_over_error', 'phot_g_mean_mag', 'phot_bp_n_obs', 'phot_bp_mean_flux', 'phot_bp_mean_flux_error', 'phot_bp_mean_flux_over_error', 'phot_bp_mean_mag', 'phot_rp_n_obs', 'phot_rp_mean_flux', 'phot_rp_mean_flux_error', 'phot_rp_mean_flux_over_error', 'phot_rp_mean_mag', 'phot_bp_rp_excess_factor', 'phot_proc_mode', 'bp_rp', 'bp_g', 'g_rp', 'radial_velocity', 'radial_velocity_error', 'rv_nb_transits', 'rv_template_teff', 'rv_template_logg', 'rv_template_fe_h', 'l', 'b', 'ecl_lon', 'ecl_lat', 'priam_flags', 'teff_val', 'teff_percentile_lower', 'teff_percentile_upper', 'a_g_val', 'a_g_percentile_lower', 'a_g_percentile_upper', 'e_bp_min_rp_val', 'e_bp_min_rp_percentile_lower', 'e_bp_min_rp_percentile_upper', 'flame_flags', 'radius_val', 'radius_percentile_lower', 'radius_percentile_upper', 'lum_val', 'lum_percentile_lower', 'lum_percentile_upper', 'gaia_astrometric_params', 'obj_name', 'obj_id', 'ra_2', 'dec_2', 'ra_error_2', 'dec_error_2', 'epoch_mean', 'zone_id', 'obj_info_flag', 'quality_flag', 'designation', 'phot_variable_flag', 'datalink_url', 'epoch_photometry_url', 'original_ext_source_id'])
	
	if writetext:
	    if os.path.exists(catname+'.txt'):
	        os.remove(catname+'.txt')
	    ascii.write(p, catname+'.txt')
	if writeldac:
		if os.path.exists(catname + '.ldac'):
			os.remove(catname + '.ldac')
		save_table_as_ldac(p, catname + '.ldac')
   
	

list_files = glob.glob('*.fits')
image=fits.open(list_files[0])
data=image[0].data
header=image[0].header
w=WCS(header)
[RA,DEC]=w.all_pix2world(data.shape[0]/2, data.shape[1]/2,1)
print(RA)
print(DEC)

#c=SkyCoord(ra=RA,dec=DEC,unit=(u.hourangle, u.deg))

#ra= c.ra.deg
#dec=c.dec.deg
#print(ra)
#print(dec)

make_gaia_catalog(RA, DEC, 0.1, 10, 20, catname='gaiacatalog')	
	


