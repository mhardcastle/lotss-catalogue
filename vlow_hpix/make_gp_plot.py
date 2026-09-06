import matplotlib.pyplot as plt
import aplpy

from matplotlib import rc

fontsize=9
rc('font',family='serif',serif=['Times'],size=fontsize)
rc('text', usetex=True)

gc=aplpy.FITSFigure('feathered-4096-reproj.fits')
gc.recenter(53,0,width=54,height=12)
gc.show_colorscale(cmap='inferno',stretch='arcsinh',vmin=1000)
gc.show_regions('greenplus.reg')
gc.add_colorbar()
gc.colorbar.set_location('top')
gc.colorbar.set_pad(0)
gc.colorbar.set_axis_label_text('144-MHz brightness temperature (K)')
gc.save('gc.pdf')
plt.show()
