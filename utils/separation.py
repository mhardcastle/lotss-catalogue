import numpy as np

def separation(c_ra,c_dec,ra,dec):
    # all values in degrees
    return np.sqrt((np.cos(c_dec*np.pi/180.0)*(ra-c_ra))**2.0+(dec-c_dec)**2.0)
