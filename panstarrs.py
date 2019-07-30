from astroquery.vizier import Vizier
import astropy.units as u
import astropy.coordinates as coord
 
Ra = input("Enter RightAscension in degress:")
Dec = input("Enter Declination in degrees:")
Rad = input("Radius of Search in degrees:")

# defining a function for querying panstarrs catalog

def panstarrs_query(ra_deg, dec_deg, rad_deg, maxmag=20, maxsources=10000):

    vquery = Vizier(columns=['*'],column_filters={"gmag": ("<%f" % maxmag)},row_limit=maxsources)
                             
    
    field = coord.SkyCoord(ra=ra_deg, dec=dec_deg, unit=(u.deg, u.deg), frame='icrs')
    return vquery.query_region(field,width=("%fd" % rad_deg), catalog="II/349/ps1")[0]
    
print(panstarrs_query(Ra, Dec, Rad))
