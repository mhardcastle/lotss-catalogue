import matplotlib.pyplot as plt
import aplpy

from matplotlib import rc

fontsize=18
rc('font',family='serif',serif=['Times'],size=fontsize)
rc('text', usetex=True)

fig=plt.figure(figsize=(15,15))
f1=aplpy.FITSFigure('nps-1.fits',figure=fig,subplot=[0.1,0.1,0.32,0.8])
f1.show_colorscale(cmap='inferno',stretch='arcsinh',vmin=750,vmax=1500)
f2=aplpy.FITSFigure('nps-0.fits',figure=fig,subplot=[0.50,0.1,0.38,0.8])
f2.show_colorscale(cmap='inferno',stretch='arcsinh',vmin=750,vmax=1500)
f2.add_colorbar()
f2.colorbar.set_location('right')
f2.colorbar.set_pad(0)
f2.colorbar.set_axis_label_text('144-MHz brightness temperature (K)')
plt.savefig('nps.pdf')

