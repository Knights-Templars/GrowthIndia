import argparse
import panquery as pq
from astropy.table import Table

parser = argparse.ArgumentParser()
parser.add_argument('filename',help="RA DEC list file, should be csv with tileid, ra, dec as three columns")
args = parser.parse_args()

coordlist = Table.read(args.filename)

for row in coordlist:
	basename = 'ps_'+str(row['tile'])
	ra = row['ra']
	dec = row['dec']
	while(1):
		try:
			pq.savefits(str(ra),str(dec),basename,size=6000)
			break
		except:
			print "Connection error...retrying..."

