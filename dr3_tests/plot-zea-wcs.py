import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from astropy.wcs import WCS
from reproject import reproject_from_healpix
from astropy import units as u
import os, sys
from astropy.io import fits

plt.rcParams['font.family'] = 'DejaVu Sans'

plt.rcParams.update({
    'font.size': 12,
    'axes.labelsize': 12,
    'xtick.labelsize': 11,
    'ytick.labelsize': 11,
    'axes.titlesize': 13,
})

# -------------------- INPUT MAP --------------------
NSIDE = 512
NPIX = hp.nside2npix(NSIDE)
print('NPIX =', NPIX)

m = np.full(NPIX, hp.UNSEEN, dtype=float)

outtextfile = f'healpix-{NPIX}-allpixels-better.txt'
if not os.path.exists(outtextfile):
    print('No existing healpix file found. Exiting.')
    sys.exit(0)

with open(outtextfile) as fh:
    for line in fh:
        pixid, sensit, hpra, hpdec = map(float, line.strip().split(','))
        m[int(pixid)] = sensit

m[np.where(m==0.0)]= 50.0  # set zero values to white

# -------------------- ZEA WCS --------------------
nx = ny = 4000
RA0  = 180.0            # central RA (deg) – change if you like
DEC0 = 90.0             # north pole at centre

pixscale = 90.0 / (min(nx, ny) / 2.0)   # deg/pixel (0.045 for 4000x4000)

w = WCS(naxis=2)
w.wcs.ctype = ['RA---ZEA', 'DEC--ZEA']
w.wcs.crval = [RA0, DEC0]
w.wcs.crpix = [nx/2.0, ny/2.0]
w.wcs.cdelt = [+pixscale, -pixscale]    # RA clockwise, Dec upwards

# -------------------- REPROJECT --------------------
array, footprint = reproject_from_healpix(
    (m, 'ICRS'),            # change to 'GALACTIC' if needed
    w.to_header(),
    shape_out=(ny, nx),
    nested=False,
    order='bilinear'
)

# -------------------- MASK TO Dec >= 0° --------------------
yy, xx = np.mgrid[0:ny, 0:nx]
ra_deg, dec_deg = w.pixel_to_world_values(xx, yy)
mask = np.isfinite(array) & (footprint > 0) & (dec_deg >= 0.0)
arr = np.ma.masked_where(~mask, array)

# -------------------- PLOT --------------------
# colormap: masked -> transparent; values < vmin -> white
cmap = plt.cm.get_cmap('cubehelix').copy()
#cmap.set_bad((1, 1, 1, 0))   # masked pixels transparent
#cmap.set_under((1, 1, 1, 1)) # under-range pixels white

vmin, vmax = -0.1, 1.5

fig = plt.figure(figsize=(10, 10), facecolor='white')
ax = fig.add_subplot(111, projection=w, facecolor='white')

im = ax.imshow(arr, origin='lower', cmap=cmap, vmin=vmin, vmax=vmax)

# crop exactly to Dec >= 0° region
ys, xs = np.where(mask)
ax.set_xlim(xs.min()-0.5, xs.max()+0.5)
ax.set_ylim(ys.min()-0.5, ys.max()+0.5)

# -------------------- RA/DEC AXES --------------------
ax.coords['ra'].set_format_unit(u.deg)
ax.coords['ra'].set_ticks(spacing=20*u.deg)
ax.coords['ra'].set_major_formatter('d')
ax.coords['ra'].set_ticklabel_position('b')
ax.coords['ra'].set_axislabel_position('b')
ax.coords['ra'].set_axislabel('RA (deg)')

dec_vals = (np.arange(0, 91, 10) * u.deg)
ax.coords['dec'].set_ticks(values=dec_vals)
ax.coords['dec'].set_major_formatter('d')
ax.coords['dec'].set_ticklabel_visible(False)
ax.coords['dec'].set_axislabel_position('l')
ax.coords['dec'].set_axislabel('Dec (deg)')  
# grid
ax.grid(color='0.7', alpha=0.5, lw=0.8)

RA_label = 270.0
cx, cy = nx/2.0, ny/2.0
inward_frac = 0.01  # 

for d in range(0, 91, 10):
    # world -> pixel for the point on the circle
    px, py = w.world_to_pixel_values(RA_label, d)
    # nudge inward toward the pole (centre of the projection)
    px2 = px + (cx - px) * inward_frac
    py2 = py + (cy - py) * inward_frac
    ax.text(
        px2, py2, f"{d}°",
        transform=ax.get_transform('pixel'),
        ha='left', va='center', rotation=0,
        fontsize=11, color='0.25',
        bbox=dict(facecolor='white', edgecolor='none', alpha=0.85, pad=0.2)
    )

ax.set_xlabel('RA (deg)')
ax.set_ylabel('Dec (deg)') 

cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
cbar.set_label('mJy/beam')

# balanced margins — left and right about the same
plt.subplots_adjust(left=0.07, right=0.92, bottom=0.08, top=0.93)

# put "Dec (deg)" a bit closer to the plot, vertically centered
fig.text(0.05, 0.5, 'Dec (deg)',
         rotation=90, va='center', ha='center', fontsize=11)



# save
plt.savefig('zea-wcs-projection-clean.png', dpi=500, facecolor='white', bbox_inches='tight')
print('Saved: zea-wcs-projection-clean.png')
