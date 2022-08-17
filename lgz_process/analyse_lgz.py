from astropy.table import Table
import numpy as np

t=Table.read('LGZ-cat.fits')

print 'Total number of rows:    ',len(t)
print 'Total with opt id:       ',np.sum((t['OptID_Name']!="None") & (t['OptID_Name']!="Mult"))
print 'Total blend>0.5:         ',np.sum((t['Blend_prob']>0.5))
print 'Total zoom>0.5:          ',np.sum((t['Zoom_prob']>0.5))
print 'Total with >1 bad clicks:',np.sum((t['Badclick']>1))
print 'Total with >1 bad + NoID:',np.sum((t['Badclick']>1) & (t['OptID_Name']=='None'))
