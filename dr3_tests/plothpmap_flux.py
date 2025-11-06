import healpy as hp
import matplotlib.pyplot as plt
import sys

hmap=hp.fitsfunc.read_map(sys.argv[1])

from healpy.newvisufunc import projview
projview(
    hmap,
    coord=["E"],
    graticule=True,
    graticule_labels=True,
    unit="RMS",
    xlabel="Right ascension (ICRS)",
    ylabel="Declination (ICRS)",
    cb_orientation="vertical",
    min=0,
    max=1.5e-3,
    latitude_grid_spacing=30,
    projection_type="aitoff",
    rot=[180,0,0],custom_xtick_labels=["300°","240°","180°","120°","60°"],
    cmap='viridis'
)
plt.tight_layout()
if len(sys.argv)>2:
    plt.savefig(sys.argv[1].replace('.fits','.pdf'))
else:
    plt.show()
