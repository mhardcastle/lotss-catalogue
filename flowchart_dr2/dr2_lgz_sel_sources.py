from astropy.table import Table
import numpy as np    

t = Table.read('LoTSS_DR2_rolling.srl_0h.sorted_step3.fits')


'''
this script doesn't actually do anything: just notes for the record:
- the selection for the 0hr field for the next set of lgz sources is the following FC_flag3 == 3,7,13
'''

t_flag7 = t[t['FC_flag3']==7]

t_flag3 = t[t['FC_flag3']==3]


t_flag13 = t[t['FC_flag3']==13]
