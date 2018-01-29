# Process the blend output files

from astropy.table import Table,vstack
import numpy as np
import sys
import os
from process_lgz import sourcename,Make_Shape

if __name__=='__main__':

    tname=sys.argv[1]
    lname=sys.argv[1].replace('.fits','-list.txt')
    t=Table.read(tname)
    mask=(t['ID_flag']<61) | (t['ID_flag']>63)
    tout=t[mask]
    tb=t[~mask]
    ct=Table.read('LOFAR_HBA_T1_DR1_merge_ID_v0.9.comp.fits')
    gt=Table.read('lofar_gaus_pw.fixed.fits')

    for i,r in enumerate(tb):
        name=r['Source_Name']
        gc=0
        gl=[]
        print name
        ctf=(ct['New_Source_Name']==name)
        ctfl=ct[ctf]
        print '... has',len(ctfl),'components'
        for j,c in enumerate(ctfl):
            print '    Component',j,'has type',c['S_Code']
            gtf=(gt['Source_Name']==c['Source_Name'])
            gtfl=gt[gtf]
            print '    ... and contains',len(gtfl),'Gaussians'
            gl.append(gtfl)
        rgt=vstack(gl)
        print '... altogether',len(rgt),'Gaussians'
        rgt['Source']=0
        dc=rgt['Maj']**2.0-6.0**2.0
        dc[dc<0]=0
        rgt['DC_Maj']=np.sqrt(dc)
        dc=rgt['Min']**2.0-6.0**2.0
        dc[dc<0]=0
        rgt['DC_Min']=np.sqrt(dc)
        
        filename=name+'.txt'
        lines=[l.rstrip() for l in open(filename).readlines()]

        # parse the file
        if lines[0]=='## Flagged':
            r['ID_flag']=610
            tout.add_row(r)
        elif lines[0]=='## Unchanged':
            r['ID_flag']=600
            tout.add_row(r)
        elif lines[0]=='## Components':
            print 'Output file to be processed!'
            for i in range(len(rgt)):
                bits=lines[i+1].split()
                rgt['Source'][i]=int(bits[1])

            optid={}
            for l in lines[len(rgt)+3:]:
                bits=l.split()
                id=int(bits[0])
                ra=float(bits[1])
                dec=float(bits[2])
                optid[id]=(ra,dec)
                
            if np.all(rgt['Source']==1):
                print 'Components unchanged'
                # in this case only the opt ID has changed
                r['ID_flag']=601
                if 1 in optid:
                    ra,dec=optid[1]
                    r['ID_ra']=ra
                    r['ID_dec']=dec
                    r['ID_name']='Altered'
                else:
                    # opt id was removed
                    r['ID_name']=""
                    r['ID_ra']=np.nan
                    r['ID_dec']=np.nan
                tout.add_row(r)
            else:
                print "It's complicated"
                # This means that the source has been split into more
                # than one component, possibly each with an optical
                # ID. We need to create one new line in the output
                # table for each input line
                ss=set(rgt['Source'])
                for s in ss:
                    if s==0: continue
                    clist=rgt[rgt['Source']==s]
                    print 'Doing source',s,'with',len(clist),'Gaussians'
                    r['ML_LR']=np.nan
                    if len(clist)==1:
                        # New source consists of only one Gaussian. Build the table entry from that
                        rg=clist[0]
                        copy=['RA','E_RA','DEC','E_DEC','Peak_flux','E_Peak_flux','Total_flux','E_Total_flux','Maj','E_Maj','Min','E_Min','DC_Maj','DC_Min','PA','E_PA','Isl_rms','Mosaic_ID']
                        for k in copy:
                            r[k]=rg[k]
                        r['Source_Name']=sourcename(r['RA'],r['DEC'])
                        r['ID_flag']=602
                        r['S_Code']='S'
                        r['ML_LR']=np.nan
                        r['LGZ_Size']=np.nan
                        r['LGZ_Width']=np.nan
                        r['LGZ_PA']=np.nan
                        r['LGZ_Assoc']=np.nan
                        r['LGZ_Assoc_Qual']=np.nan
                        r['LGZ_ID_Qual']=np.nan
                        r['DC_PA']=r['PA']
                        for k in ['Min','Maj','PA']:
                            r['E_DC_'+k]=r['E_'+k]
                    else:
                        # New source consists of several components. Code from process_lgz
                        print "It's really complicated"
                        tfluxsum=np.sum(clist['Total_flux'])
                        ra=np.sum(clist['RA']*clist['Total_flux'])/tfluxsum
                        dec=np.sum(clist['DEC']*clist['Total_flux'])/tfluxsum
                        sname=sourcename(ra,dec)
                        print '      New sourcename is',sname
                        r['ID_flag']=603
                        r['RA']=ra
                        r['DEC']=dec
                        r['Source_Name']=sname
                        r['E_RA']=np.sqrt(np.mean(clist['E_RA']**2.0))
                        r['E_DEC']=np.sqrt(np.mean(clist['E_DEC']**2.0))
                        r['Total_flux']=np.sum(clist['Total_flux'])
                        r['E_Total_flux']=np.sqrt(np.sum(clist['E_Total_flux']**2.0))
                        maxpk=np.argmax(clist['Peak_flux'])
                        r['Peak_flux']=clist[maxpk]['Peak_flux']
                        r['E_Peak_flux']=clist[maxpk]['E_Peak_flux']
                        r['S_Code']='M'
                        r['Isl_rms']=np.mean(clist['Isl_rms'])
                        r['Mosaic_ID']=clist[maxpk]['Mosaic_ID']
                        ms=Make_Shape(clist)
                        r['LGZ_Size']=ms.length()
                        r['LGZ_Width']=ms.width()
                        r['LGZ_PA']=ms.pa()
                        r['LGZ_Assoc']=len(clist)
                        r['LGZ_Assoc_Qual']=1
                        r['LGZ_ID_Qual']=1
                        r['Number_Masked']=0
                        r['Number_Pointings']=0
                        r['Masked_Fraction']=0
                        for k in ['Maj','E_Maj','Min','E_Min','DC_Maj','DC_Min','PA','E_PA']:
                            r[k]=np.nan
                    if s in optid:
                        ra,dec=optid[s]
                        r['ID_ra']=ra
                        r['ID_dec']=dec
                        r['ID_name']='Altered'
                    else:
                        r['ID_name']=""
                        r['ID_ra']=np.nan
                        r['ID_dec']=np.nan
                    tout.add_row(r)
        else:
            raise RuntimeError('Cannot parse input file...'+lines[0])

    tout.write('output.fits',overwrite=True)
    
            
        

    
