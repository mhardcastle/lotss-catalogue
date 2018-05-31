# Process the blend output files

from astropy.table import Table,vstack
import numpy as np
import sys
import os
from process_lgz import sourcename,Make_Shape
from copy import deepcopy

def add_artefact(at,r):
    new_row=deepcopy(at[0])
    for c in at.colnames:
        new_row[c]=r[c]
    new_row['Source_Name']=r['Component_Name']
    print '    *** Adding %s to artefact list ***' % new_row['Source_Name']
    at.add_row(new_row)

if __name__=='__main__':

    t=Table.read('LOFAR_HBA_T1_DR1_merge_ID_v1.2.fits')
    t['Deblended_from']='                      '
    mask=(t['ID_flag']<41) | (t['ID_flag']>42)
    tout=t[mask]
    tb=t[~mask]
    # component table. This is modified too
    ct=Table.read('LOFAR_HBA_T1_DR1_merge_ID_v1.2.comp.fits')
    ct['Deblended_from']='                      '
    mask=(ct['ID_flag']<41) | (ct['ID_flag']>42)
    ctout=ct[mask]
    # artefact table. Append dropped sources to this
    at=Table.read('LOFAR_HBA_T1_DR1_merge_ID_v1.2.art.fits')
    delcols=[]
    for c in at.colnames:
        if c not in ct.colnames:
            delcols.append(c)
    at.remove_columns(delcols)
    # Gaussian table, for components
    gt=Table.read('lofar_gaus_pw.fixed.fits')

    for i,r in enumerate(tb):
        name=r['Source_Name']
        gc=0
        gl=[]
        print name
        ctf=(ct['Source_Name']==name)
        ctfl=ct[ctf]
        print '... has',len(ctfl),'components'
        for j,c in enumerate(ctfl):
            print '    Component',j,'has type',c['S_Code']
            gtf=(gt['Source_Name']==c['Component_Name'])
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
            #tout.add_row(r) # drop these
            for cr in ctfl:
                #cr['ID_flag']=610
                add_artefact(at,cr)
        elif lines[0]=='## Unchanged':
            tout.add_row(r)
            for cr in ctfl:
                ctout.add_row(cr)
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
            if np.all(rgt['Source']==0):
                print 'No unflagged components!'
                # treat as flagged
                for cr in ctfl:
                    add_artefact(at,cr)
            elif np.all(rgt['Source']==1):
                print 'Components unchanged'
                # in this case only the opt ID has changed
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
                for cr in ctfl:
                    cr['ID_flag']=r['ID_flag']
                    ctout.add_row(cr)
            else:
                print "It's complicated"
                # This means that the source has been split into more
                # than one component, possibly each with an optical
                # ID. We need to create one new line in the output
                # table for each input line
                ss=set(rgt['Source'])
                for s in ss:
                    if s==0:
                        # implies there are Gaussians here that don't belong
                        # in any source.
                        clist=rgt[rgt['Source']==0]
                        print 'Checking source',s,'for artefact components'
                        # if there are whole components none of which
                        # are in the source, then we need to add to the
                        # artefact list
                        print clist
                        components=set(clist['Source_Name'])
                        for c in components:
                            flagged=np.sum(clist['Source_Name']==c)
                            total=np.sum(gt['Source_Name']==c)
                            if flagged==total:
                                crow=ctfl[ctfl['Component_Name']==c][0]
                                print 'Flagging component',crow['Component_Name']
                                add_artefact(at,crow)
                    clist=rgt[rgt['Source']==s]
                    print 'Doing source',s,'with',len(clist),'Gaussians'
                    r['ML_LR']=np.nan
                    r['Deblended_from']=name
                    if len(clist)==1:
                        # New source consists of only one Gaussian. Build the table entry from that
                        rg=deepcopy(clist[0])
                        copy=['RA','E_RA','DEC','E_DEC','Peak_flux','E_Peak_flux','Total_flux','E_Total_flux','Maj','E_Maj','Min','E_Min','DC_Maj','DC_Min','PA','E_PA','Isl_rms','Mosaic_ID']
                        for k in copy:
                            r[k]=rg[k]
                        sname=sourcename(r['RA'],r['DEC'])
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
                        r['RA']=ra
                        r['DEC']=dec
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
                        for k in ['Maj','E_Maj','Min','E_Min','DC_Maj','DC_Min','PA','E_PA','DC_PA','E_DC_Maj','E_DC_Min', 'E_DC_PA']:
                            r[k]=np.nan
                    r['Source_Name']=sname
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
                    # now sort out the components by converting
                    # Gaussians to new components to go in ctout
                    for g in clist:
                        c=deepcopy(ct[0])
                        copy=['RA','E_RA','DEC','E_DEC','Peak_flux','E_Peak_flux','Total_flux','E_Total_flux','Maj','E_Maj','Min','E_Min','DC_Maj','DC_Min','PA','E_PA','Isl_rms','Mosaic_ID']
                        for k in copy:
                            c[k]=g[k]
                        for k in ['Min','Maj','PA']:
                            c['E_DC_'+k]=c['E_'+k]
                        # Gaussians don't have names
                        c['Deblended_from']=g['Source_Name'] # the original component name
                        c['Component_Name']=sourcename(g['RA'],g['DEC'])
                        c['Source_Name']=sname
                        c['ID_flag']=r['ID_flag']
                        # now fix up a few other columns
                        #c['Artefact_flag']=False
                        #c['LGZ_flag']=0
                        #c['FC_flag']=0
                        c['Ng']=1
                        c['S_Code']='S'
                        #c['G_max_sep']=np.nan
                        for k,ty in ct.dtype.descr:
                            if k[0:2]=='NN' or 'LR' in k:
                                if ty=='>f8':
                                    c[k]=np.nan
                                elif ty=='>i8':
                                    c[k]=0
                                elif 'S' in ty:
                                    c[k]=''
                        ctout.add_row(c)
                    
        else:
            raise RuntimeError('Cannot parse input file...'+lines[0])

    tout.write('merge_out.fits',overwrite=True)
    ctout.sort('RA')
    ctout.write('merge_comp_out.fits',overwrite=True)
    at.sort('RA')
    at.write('merge_art_out.fits',overwrite=True)
    
    
