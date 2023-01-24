from astropy.table import Table, vstack


t0 = Table.read('LoTSS_DR2_v110.srl_0h.lr-full.sorted_step3_flux4.fits')
t13 = Table.read('LoTSS_DR2_v110.srl_13h.lr-full.sorted_step3_flux4.fits')

t = vstack([t0,t13])

t.write('LoTSS_DR2_v110.srl.lr-full.sorted_step3_flux4.fits')
