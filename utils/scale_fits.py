"""
Scale data in fits file by a factor
TODO: Inplace mode is not working
Contributed by Pepe
"""
import argparse
import os
from pathlib import Path
from astropy.io import fits
import numpy as np


def scale_fits(fits_name, scale, output=None, del_naxis=False):
    """Scale fits image.
    """
    if output is None:
        mode = 'update'
    else:
        mode = 'readonly'
    with fits.open(fits_name, mode=mode) as hdu:
        hdu[0].data = hdu[0].data / scale
        hdr = hdu[0].header
        hdr['history'] = "Data scaled by {}".format(scale)
        if len(hdu[0].data) == 4 and del_naxis:
            del hdu[0].header['NAXIS3']
            del hdu[0].header['NAXIS4']
        if output is None:
            hdu.flush()
        else:
            hdu.writeto(output,overwrite=True)

def main(args):
    for fits_file in args.fits_images:
        path = Path(fits_file.name)
        if not args.inplace:
            output = path.parent / Path(path.stem + "-scaled" + path.suffix)
            scale_fits(fits_file, args.scale, output=output)
        else:
            scale_fits(fits_file, args.scale)
        

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Scale fits images.',
        epilog='It creates a new file with the suffix "-scaled" by default.')
    parser.add_argument('fits_images', 
                        metavar='FITS',
                        nargs='+',
                        type=argparse.FileType('rb'),
                        help='Fits images to scale')
    parser.add_argument('-s', '--scale',
                        type=float,
                        required=True,
                        help='Scaling factor (the values are divided by it)')
    parser.add_argument('-i', '--inplace',
                        action='store_true',
                        help='Overwrite the input file')
    args = parser.parse_args()
    main(args)

