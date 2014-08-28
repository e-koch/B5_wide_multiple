import os
from astropy.io import fits

file_scuba2_raw='B5_450um_ext_v2.fits'
file_scuba2_out='B5_450um_ext_v2_beam.fits'
hdu = fits.open(file_scuba2_raw)
hdr =hdu[0].header
data=hdu[0].data
hdu.close()
hdr.append(('BMAJ', 9.8/3600.))
hdr.append(('BMIN', 9.8/3600.))
hdr.append(('BPA',  0.0))
os.system('rm -r '+file_scuba2_out)
fits.writeto(file_scuba2_out, data, hdr)
