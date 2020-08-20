import sys
from astropy.table import Table

filename = sys.argv[1]

tt = Table.read(filename)
if 'hdf5' in filename:
    tt.write(filename.replace('hdf5','fits'), overwrite=True)
elif 'fits' in filename:
    tt.write(filename.replace('fits','hdf5'), serialize_meta=True, overwrite=True)
