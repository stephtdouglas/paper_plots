# Module of functions to create plots for AAS 221 poster
# Created by Stephanie Douglas, 1 December 2012
# Copied from praesepe_comp.py 19 December 2012
# Copied from praesepe_comp_v2.py 21 December 2012
# Copied from aas221.py 27 February 2013
# Copied from completeness_Feb2013.py 2 March 2013
# ******************************************************************************

import numpy as np
import matplotlib.pyplot as plt
import asciitable as at
import pyfits
import get_data

kh = at.read('/home/stephanie/ptf/models/kraushillenbrand5.dat')
ksub = np.array([-23,-19,-15,-12,-11,-9,-7,-6])#arange(-15,-1,1)
kh_rpmK = (kh['Mr'] - 0.035*(kh['Mr']-kh['Mi']) - 0.007 - kh['MK'])[ksub]
kh_rpmK[0] += 0.1
kh_spt = kh['SpT'][ksub]
klen = len(kh_spt)



def getspt(spt_str):
    """
    Converts a list of spectral types into numbers for plotting

    Inputs:
        spt_str - an array of strings giving spectral types

    Outputs:
        stype_num - an array of numbers corresponding to spectral types
            ('G0' = 0)
    """

    snum = len(spt_str)

    stype = np.empty(snum,'S2')
    stype[:] = '  '
    stype_num = np.ones(snum,int)*-99

    for i in range(snum):
        stype[i] = spt_str[i].split('.')[0]
        if stype[i][0]=='G':
            stype_num[i] = int(stype[i][1])
        elif stype[i][0]=='K':
            stype_num[i] = int(stype[i][1])+10
        elif stype[i][0]=='M':
            stype_num[i] = int(stype[i][1])+20
        elif stype[i][0]=='L':
            stype_num[i] = int(stype[i][1])+30
        elif stype[i][0]=='F':
            stype_num[i] = int(stype[i][1])-10
        elif stype[i][0]=='A':
            stype_num[i] = int(stype[i][1])-20
        elif stype[i][0]=='B':
            stype_num[i] = int(stype[i][1])-30
        elif stype[i][0]=='O':
            stype_num[i] = int(stype[i][1])-40

    return stype_num

def plot_stars(data,stars,ax,cluster='P',ptype='cmd',clr='rK',pmem=False,
    ploterrs=True,obs=True):
    """
    inputs:
        data - data structure from a fits file
        stars - array defining the indices in data for the stars to plot
        cluster='P' - defaults to 'P' for Praesepe, 'H' is for Hyades
        ptype='cmd' - defaults to a color-magnitude diagram of (r-K,r)
            other options are 'Prot' which plots (color,Prot)
        clr='rK' - defaults to 'rK' for (r-K), 
            other options are 'JK' for (J-K)
    """
    for i in stars:
#        print i
        # Get the r and K magnitudes
        r = data.field('RPRIME')[i]
        rerr = data.field('RPRIME_ERR')[i]
        K = data.field('TWOMASS_K')[i]
        Kerr = data.field('TWOMASS_KERR')[i]
        J = data.field('TWOMASS_J')[i]
        Jerr = data.field('TWOMASS_JERR')[i]
        if rerr>0.5 or Kerr>0.5 or Jerr>0.5:
            print r, rerr, K, Kerr, J, Jerr
#        print r, J, K
        if cluster=='P' and r<8:
            print r,data.field('RMAG_FLAG')[i],(r-K)
        
        # The color of the marker indicates the source of the period
        # I'm going to make the markers a little bigger for rotators
#        print cluster
        if cluster=='P':
            if (data.field('PERIOD_FLAG')[i]=='P'): #PTF period
                #print 'PTF'
                my_mc = 'Red'
                mk_sz = 9
                prot = data.field('PTF_PERIOD')[i]
            elif (data.field('PERIOD_FLAG')[i]=='D'): 
                #print 'Del'
                #Delorme period
                my_mc = '#00CCFF'
                mk_sz = 9
                prot = data.field('SWASP_PERIOD')[i]
            elif (data.field('PERIOD_FLAG')[i]=='S'):
                #print 'Scholz'
                my_mc = 'DarkSlateBlue'
                mk_sz = 9
                prot  = data.field('SCHOLZ_PERIOD')[i]
            else:
                #print 'No'
                my_mc = 'Black'
                mk_sz = 8
                prot = -99.
        if cluster=='H':
            if (data.field('PERIOD_FLAG')[i]=='D'):
                #Delorme period
                my_mc = '#00CCFF'
                mk_sz = 9
                prot = data.field('DELORME_LITP')[i]
            elif (data.field('PERIOD_FLAG')[i]=='K'):
                my_mc = 'DarkGreen'
                mk_sz = 9
                prot = data.field('KUNDERT_PROT')[i]
                #print prot
            elif (data.field('PERIOD_FLAG')[i]=='R'):
                my_mc = 'Goldenrod'
                mk_sz = 9
                prot = data.field('PERIOD')[i]
                #print prot
            else:
                my_mc = 'Black'
                mk_sz = 8
                prot = -99.
        if obs:
            my_mfc = my_mc
        else:
            my_mfc = 'None'


        # If membership probability is included, 
        # then only plot members with >75% membership probability
        #print pmem
        if pmem and cluster=='P':
            minpmem=75.0
            pmem = data.field('ADAMPMEM')[i]
            if pmem<0:
                pmem=76.0
        else:
            minpmem=0.0
            pmem=100.0



#        print ptype        
        if ptype=='cmd':
            # Shape of the marker indicates the source of the spectrum
            if obs==False:
                my_sh = 'o'
                mk_sz = 5
            elif cluster=='H':
                my_sh = '*'
                mk_sz += 3
            elif data.field('MDM_SPECMATCH')[i]>0:
                my_sh = '*'
                mk_sz += 1
            elif data.field('MAGE_SPECMATCH')[i]>0:
                my_sh = 'p'
            elif data.field('HYDRA_SPECMATCH')[i]>0:
                my_sh = 'D'
                mk_sz=mk_sz-1
            elif data.field('KAFKAMATCH')[i]>0:
                my_sh = 'v'
            elif data.field('ALLEN_SPECMATCH')[i]>0:
                my_sh = '^'
            elif data.field('SDSS_SPECMATCH')[i]>0:
                my_sh = 's'
            #print 'cmd',my_sh
        elif ptype=='Prot':
            if cluster=='H':
                my_sh = 'd'
            else:
                my_sh = '*'
                mk_sz += 1
        else:
            my_sh = '*'
            mk_sz += 1
        
#        print 'got type'

        #print pmem,minpmem
        if (ptype=='cmd' and (
            (pmem>minpmem) or (cluster=='H'))):
            #print pmem,my_sh
            if clr=='rK':
                xcolor = r-K
                ymag = r
            elif clr=='JK':
                xcolor = J-K
                ymag = J
            if (clr=='rK') and ploterrs:
                xcolor_err = np.sqrt(rerr**2+Kerr**2)
                ymag_err = rerr
            elif (clr=='JK') and ploterrs:
                xcolor_err = np.sqrt(Jerr**2+Kerr**2)
                ymag_err = Jerr
            else:
                xcolor_err=0
                ymag_err=0
            #print xcolor, ymag, xcolor_err, ymag_err
            #print my_mc, my_mfc, my_sh, mk_sz, my_al
            ax.errorbar(xcolor, ymag, yerr=ymag_err, xerr=xcolor_err,
                ecolor='Black',capsize=0, 
                mec=my_mc, 
                mfc=my_mfc, 
                marker=my_sh, mew=0.75,
                ms=mk_sz)#, alpha=my_al)
        elif ptype=='Prot' and prot>0:
#            print ptype,clr
            if clr=='rK':
                xcolor = (r-K)
                xcolor_err = np.sqrt(rerr**2+Kerr**2)
            elif clr=='JK':
                xcolor = (J-K)
                xcolor_err = np.sqrt(Jerr**2+Kerr**2)
            if ploterrs==False:
                xcolor_err=0
#            print xcolor, xcolor_err
            ax.errorbar(xcolor, prot, xerr=xcolor_err,
                ecolor='Black',capsize=0,
                mec=my_mc, mfc=my_mfc, marker=my_sh, ms=mk_sz)#,alpha=my_al)
        elif (ptype!='Prot') and (ptype!='cmd'):
            print 'plot type %s not understood' % str(ptype)
        #print 'plotted!'

def plot_notobs_wrap(data,obs,obs_nonrot,obs_rot,cluster,ax1,
    p_type,colr='rK',p_mem=False,ploterrs=True,outfile=None,
    newfig=True):

    not_obs = np.ones(len(data))
    not_obs[obs] = 0
    not_obs = np.where(not_obs)[0]

    rot = np.where(data.field('ANYPERIOD')==1)[0]
    nonrot = np.where(data.field('ANYPERIOD')==0)[0]

    no_nr = np.intersect1d(not_obs,nonrot)
    no_rot = np.intersect1d(not_obs,rot)

    if newfig:
        plt.figure(figsize=(9,7))
        ax1 = plt.subplot(111)
    if p_type=='cmd':
        plot_stars(data,no_nr,ax1,cluster,p_type,colr,p_mem,ploterrs,obs=False)
#        print 'Plotted non-rotators'
    plot_stars(data,no_rot,ax1,cluster,p_type,colr,p_mem,ploterrs,obs=False)
#    print 'Plotted rotators'
    plot_stars_wrap(data,obs_nonrot,obs_rot,cluster,p_type,colr,p_mem,
        ploterrs,outfile,newfig=False,ax=ax1)


def plot_stars_wrap(data,obs_nonrot,obs_rot,cluster,
    p_type,colr='rK',p_mem=False,ploterrs=True,
    outfile=None,newfig=True,ax=None):
    """
    """

    if newfig:
        plt.figure(figsize=(9,7))
        ax = plt.subplot(111)
    ax.minorticks_on()
    ax.tick_params(which='both',width=1,labelsize='large')
#    print 'Axes Set'

    if p_type=='cmd':
        plot_stars(data,obs_nonrot,ax,cluster,p_type,colr,p_mem,ploterrs)
#        print 'Plotted non-rotators'
    plot_stars(data,obs_rot,ax,cluster,p_type,colr,p_mem,ploterrs)
#    print 'Plotted rotators'


    if p_type=='cmd':
        if colr=='rK':
            ax.set_xlim(0,7)
            
            ax.set_xlabel('(r\'-K)',fontsize='x-large')
            ax.set_ylabel('r\' ',fontsize='x-large')
            tx = 0.5
            
            if cluster=='P':
                ty = 19
                ax.set_ylim(21,6)
            else:
                ty = 16
                ax.set_ylim(18,3)                
        elif colr=='JK':
            ax.set_xlim(0,1.1)
            ax.set_xlabel('J-K',fontsize='x-large')
            ax.set_ylabel('J',fontsize='x-large')
            tx = 0.1
            if cluster=='P':
                ty = 14
                ax.set_ylim(16,7)
            else:
                ty = 11
                ax.set_ylim(13,4)
        ax.text(tx,ty-0.6,'Periods measured by',fontsize='large')
        if cluster=='P':
            ax.text(tx,ty,'Agueros+ 2011',color='Red',fontsize='large')
            ax.plot((0.955,1.025),(18.55,18.55),'rs',ms=0.75,mec='r')
            #ax.plot(tx+0.048,ty-0.23,'r.',ms=2)
            #ax.plot(tx+0.052,ty-0.23,'r.',ms=2)
            ax.text(tx,ty+0.6,'Delorme+ 2011',color='#00CCFF',fontsize='large')
            ax.text(tx,ty+1.2,'Scholz+ 2007, 2011',color='DarkSlateBlue',fontsize='large')
            symbols = ('*','D','p','v','^','s','o')
            sizes = (9,7,8,8,8,8,5)
            labels = ('ModSpec','Hydra','MagE','Kafka & Honeycutt 2006',
                'Allen & Strom 1995','SDSS','Not observed')
        elif cluster=='H':
            ax.text(tx,ty,'Kundert+ in prep',color='DarkGreen',fontsize='large')
            ax.text(tx,ty+0.6,'Delorme+ 2011',color='#00CCFF',fontsize='large')
            ax.text(tx,ty+1.2,'Radick+ 1987, 1995',color='Goldenrod',fontsize='large')
            symbols = ('*','o')
            sizes = (9,5)
            labels = ('ModSpec','Not observed')
        lnum = len(symbols)
        for l in range(lnum):
            if symbols[l]=='o':
                face='None'
            else:
                face='Black'
            ax.plot(-10,-10,mec='Black',mfc=face,marker=symbols[l],
                ms=sizes[l], mew=0.75,label=labels[l],ls='None')
        if cluster=='P':
            ax.legend(numpoints=1,handletextpad=0.3,
                handlelength=1,borderaxespad=0.2,prop={'size':12})
    elif p_type=='Prot':
        ax.semilogy()
#        print 'log'
        ax.set_ylim(0.1,40)
        ax.set_ylabel('Period (d)')
#        print 'labeled y'
        if colr=='rK':
#            print colr
            ax.set_xlim(0,7)
            ax.set_xlabel('(r\'-K)')
        elif colr=='JK':
            ax.set_xlim(0,1.2)
            ax.set_xlabel('(J-K)')
        if cluster=='P':
#            print cluster, colr
            ax.plot(-2,999,'*',mfc='Red',mec='Red',ms=9,label='PTF Rotators')
            ax.plot(-2,999,'*',mfc='#00CCFF',ms=9,
                mec='#00CCFF',label='Delorme+')
            ax.plot(-2,999,'*',mfc='DarkSlateBlue',mec='DarkSlateBlue',
                ms=8,label='Scholz+')
        elif cluster=='H':
#            print cluster, colr
            ax.plot(-2,999,'d',
                mfc='#00CCFF',mec='#00CCFF',ms=8,label='Delorme+')
            ax.plot(-2,999,'d',
                mfc='DarkGreen',mec='DarkGreen',ms=8,label='ASAS')
        ax.legend(loc=3,numpoints=1,scatterpoints=1,frameon=False)

    if outfile!=None:
        ax.savefig(outfile)



def sidebyside():#p_type,colr='rK',p_mem=False,ploterrs=True):
    """
    """
    hdat, hobs, hobs_nr, hobs_r = get_data.get_data('H')
    pdat, pobs, pobs_nr, pobs_r = get_data.get_data('P')

    plt.figure(figsize=(13,5.5))
    axh = plt.subplot(122)
    axp = plt.subplot(121)

    plot_notobs_wrap(hdat,hobs,hobs_nr,hobs_r,'H',axh,'cmd','rK',
        p_mem=False,ploterrs=False,newfig=False)
    texty = 2.75
    for i in range(klen):
        axh.text(kh_rpmK[i],texty,kh_spt[i],fontsize='large')
    plot_notobs_wrap(pdat,pobs,pobs_nr,pobs_r,'P',axp,'cmd','rK',
        p_mem=False,ploterrs=False,newfig=False)
    texty = 5.75
    for i in range(klen):
        axp.text(kh_rpmK[i],texty,kh_spt[i],fontsize='large')
    axh.tick_params(which='both',top=False)
    axp.tick_params(which='both',top=False)
    axh.text(5.75,3.75,'Hyades',fontsize='large')
    axp.text(5.4,13,'Praesepe',fontsize='large')

    plt.savefig('paper_completeness.png')
    plt.savefig('paper_completeness.ps',orientation='landscape')

sidebyside()
