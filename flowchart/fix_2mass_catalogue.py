#!/usr/bin/python

'''
fix_2mass_catalogue
apply some manual fixes to some 2MASX sources known to have bad sizes
'''

from astropy.table import Table
import os
import argparse


def fix_2mass_catalogue(filein, fileout, clobber=False):
    ## get 2MASX information
    xsc = Table.read(filein)
    
    # fix known issue with these sources - has bad size, use sdss size!
    xsc['r_ext'][xsc['designation']=='13174820+4702571    '] = 20. 
    xsc['r_ext'][xsc['designation']=='11244075+5218078    '] = 20. 
    xsc['r_ext'][xsc['designation']=='11435865+5600295    '] = 12. 
    xsc['r_ext'][xsc['designation']=='11095404+4859120    '] = 13. 
    xsc['r_ext'][xsc['designation']=='11393585+5555286    '] = 13. 


    if os.path.exists(fileout):
        if clobber:
            os.remove(fileout)
            print 'Output file, {f:s}, exists: will overwrite it'.format(f=fileout)
            xsc.write(fileout)
        else:
            print 'Output file, {f:s}, exists and clobber is set to False. Re-run with --clobber flag to overwrite it'.format(f=fileout)
    else:
        xsc.write(fileout)
    
    
    
if __name__=='__main__':
    
    parser = argparse.ArgumentParser(description='Fix 2MASX catalogue')

    parser.add_argument('-c','--clobber', dest='clobber', action='store_true', help='Overwrite existing file')
    parser.add_argument('xsc_file', help='Input 2MASX catalogue')
    parser.add_argument('xsc_file_out', help='Output 2MASX catalogue')
    args = parser.parse_args()

    
    fix_2mass_catalogue(args.xsc_file, args.xsc_file_out, clobber=args.clobber)