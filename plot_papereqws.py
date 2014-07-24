
from scipy.interpolate import interp1d
import numpy as np
import praesepe_comp as pc
import get_data 
import asciitable as at
import matplotlib.pyplot as plt
#import plot_grid as pg

kh = at.read('/home/stephanie/ptf/models/kraushillenbrand5.dat')
ksub = np.array([-15,-12,-11,-9,-7,-6])#arange(-15,-4,1)
kh_rpmK = (kh['Mr'] - 0.035*(kh['Mr']-kh['Mi']) - 0.007 - kh['MK'])[ksub]
kh_spt = kh['SpT'][ksub]
klen = len(kh_spt)


def simple(eqw,ueqw,chi):
    # 2sigma upper limits
    ulim_eqw = eqw+ueqw*2.0
    err_ulim_lha = (eqw+ueqw*2.0)*chi
    ulim_lha = err_ulim_lha*0.8
    uperr_ulim_lha = abs(err_ulim_lha - ulim_lha)
    dnerr_ulim_lha = ulim_lha*0.2#copy(uperr_ulim_lha)

    #print ulim_eqw[0:10]
    #print ulim_lha[0:10]
    #print np.where(((ulim_lha-dnerr_ulim_lha)<0) & (ulim_lha>0))[0]


    ulim_lha[(eqw>0) & ((eqw-ueqw)>0)] = -99.
    uperr_ulim_lha[(eqw>0) & ((eqw-ueqw)>0)] = -99.
    dnerr_ulim_lha[(eqw>0) & ((eqw-ueqw)>0)] = -99.

    return ulim_lha, uperr_ulim_lha, dnerr_ulim_lha


def plotit(dat,obs,ax1,plot_color,plot_marker,plot_label,spt_field,
    snum_field,offset,to_plot):
    ldat = len(dat)
    
    num_spt = pc.getspt(dat.field(spt_field))
    rmK = dat.field('RPRIME_K')
    rmKerr = dat.field('RPRIME_K_ERR')
    if len(dat)>1000:
        pmem = dat.field('ADAMPMEM')
        dkeys = ['mdm', 'kafka', 'hydra', 'allen', 'mdm2', 'sdss', 'mage']
    else:
        pmem = np.ones(len(dat))*99.0
        pmem[dat.field('ROESER_DISTANCE')<=9.0] = 99.0
        pmem[(dat.field('ROESER_DISTANCE')>9.0) & 
             (dat.field('ROESER_DISTANCE')<=18.0)] = 92.5
        pmem[dat.field('ROESER_DISTANCE')>18.0] = 70.0
        dkeys = ['mdm','mdm2']
    
    avg_eqw, unc_eqw = dat.field('AVG_EQW'),dat.field('AVG_EQW_ERR')
    avg_ll, unc_ll = dat.field('AVG_LHA'),dat.field('AVG_LHA_ERR')
    chi = dat.field('CHI')


    msize = 6
    if plot_marker=='*':
        msize += 3

    pmem_threshold = 70.0
    if to_plot=='color':
        good_stars = np.where(((pmem>=pmem_threshold) | (pmem<0))
            & (abs(unc_eqw)<10)
            & (abs(rmKerr)<1))[0]
        ax1.errorbar(rmK[good_stars],-1*avg_eqw[good_stars],
#            xerr=rmKerr[good_stars],
            yerr=unc_eqw[good_stars],mfc=plot_color,
            marker=plot_marker,mec='None',label=plot_label,lw=0,elinewidth=1,
            ms=msize,capsize=0,ecolor=plot_color)
    elif to_plot=='lha':
        good_stars = np.where(((pmem>=pmem_threshold) | (pmem<0))
            & (abs(unc_ll)<10)
            & (abs(rmKerr)<1) & ((avg_eqw-unc_eqw)>0))[0]
        ax1.errorbar(rmK[good_stars],avg_ll[good_stars],
#            xerr=rmKerr[good_stars],
            yerr=unc_ll[good_stars],
            mfc=plot_color,
            marker=plot_marker,mec='None',label=plot_label,lw=0,elinewidth=1,
            ms=msize,capsize=0,ecolor=plot_color)
        ulim_lha, uperr_ulim_lha, dnerr_ulim_lha = simple(avg_eqw,unc_eqw,chi)
        good_stars = np.where((pmem>pmem_threshold) & (abs(avg_eqw)<10)
            & (abs(rmKerr)<1) & ((avg_eqw-unc_eqw)<0) 
            & (dat.field('BINARY')==0))[0]
        pup = (np.ones(2*len(good_stars))*1e-6).reshape(2,-1)
        pup[1] = uperr_ulim_lha[good_stars]
        pup[0] = dnerr_ulim_lha[good_stars]

        if plot_color=='DarkBlue':
            error_color='#B2B2DC'#'#4D4DAE'
        else:
            error_color='#FFA280'#'#FF6A33'
        ax1.errorbar(rmK[good_stars],
            ulim_lha[good_stars],
            pup,ecolor=error_color,
            color=error_color,fmt=None,lolims=True)
#        print rmK[good_stars]
#        print ulim_lha[good_stars]
#        print pup


    elif to_plot=='lha_noerrs':
        good_stars = np.where((pmem>pmem_threshold) & (abs(unc_ll)<10)
            & (abs(rmKerr)<1))[0]
        ax1.errorbar(rmK[good_stars],avg_ll[good_stars],
#            xerr=rmKerr[good_stars],
            yerr=unc_ll[good_stars]*0.0,
            mfc=plot_color,
            marker=plot_marker,mec='None',label=plot_label,lw=0,elinewidth=1,
            ms=msize,capsize=0,ecolor=plot_color)

# r-K color
color_f = plt.figure(figsize=(9,8))
ax = plt.subplot(111)
ax.set_xlabel('(r\'-K)',fontsize='x-large')
ax.set_xlim(1.5,6.1)
ax.tick_params(which='both',width=2,labelsize='x-large',top=False)
ax.set_ylabel(r'H$\alpha$ EqW ($\AA$)',fontsize='x-large')

pdat,pobs,pobs_nr,pobs_r = get_data.get_data('P')
plotit(pdat,pobs,ax,
    'DarkBlue',#'#0099CC',
    'o','Praesepe','ADAMSPT','NUM_SPECTRA',0,'color')
hdat,hobs,hobs_nr,hobs_r = get_data.get_data('H')
plotit(hdat,hobs,ax,
    'OrangeRed',#'#FFAA00',
    'D','Hyades','MDM_SPEC_ADAMSPT','MDM_SPECMATCH',0,'color')
ax.legend(numpoints=1,handletextpad=0.2,
    handlelength=1,borderaxespad=0.2,loc=2,frameon=False)
ax.set_ylim(1.5,-15)
texty = -15.25
for i in range(klen):
    ax.text(kh_rpmK[i],texty,kh_spt[i],fontsize='large')
color_f.savefig('papereqws.png')
color_f.savefig('papereqws.eps',bbox_inches='tight')



lha_f = plt.figure(figsize=(9,8))
ax = plt.subplot(111)
ax.set_xlabel('(r\'-K)',fontsize='x-large')
ax.set_xlim(1.5,6.1)
ax.tick_params(which='both',width=2,labelsize='x-large',top=False)
ax.set_ylabel(r'$L_{H\alpha}/L_{bol}$',fontsize='xx-large')
ax.set_yscale('log')

pdat,pobs,pobs_nr,pobs_r = get_data.get_data('P')
plotit(pdat,pobs,ax,
    'DarkBlue',#'#0099CC',
    'o','Praesepe','ADAMSPT','NUM_SPECTRA',0,'lha')
hdat,hobs,hobs_nr,hobs_r = get_data.get_data('H')
plotit(hdat,hobs,ax,
    'OrangeRed',#'#FFAA00',
    'D','Hyades','MDM_SPEC_ADAMSPT','MDM_SPECMATCH',0,'lha')
ax.legend(numpoints=1,handletextpad=0.2,
    handlelength=1,borderaxespad=0.2,loc=2,frameon=False)

legendx = 1.64
ltextx = 1.73
p_y = 4.1e-4
ax.errorbar(legendx,p_y,p_y*0.15,lolims=True,ecolor='#B2B2DC',fmt=None)
ax.text(ltextx,p_y*0.9,"Praesepe upper limit",fontsize='large')
h_y = 2.8e-4
ax.errorbar(legendx,h_y,h_y*0.15,lolims=True,ecolor='#FFA280',fmt=None)
ax.text(ltextx,h_y*0.9,"Hyades upper limit",fontsize='large')

ax.set_ylim(5e-7,1e-3)
texty = 1.1e-3
for i in range(klen):
    ax.text(kh_rpmK[i],texty,kh_spt[i],fontsize='large')
lha_f.savefig('paperlhalbol.png')
lha_f.savefig('paperlhalbol.eps',bbox_inches='tight')

"""
lha_f = plt.figure(figsize=(9,8))
ax = plt.subplot(111)
ax.set_xlabel('(r\'-K)',fontsize='x-large')
ax.set_xlim(2,6.5)
ax.tick_params(which='both',width=2,labelsize='x-large')
ax.set_ylabel(r'$L_{H\alpha}/L_{bol}$',fontsize='xx-large')
ax.set_yscale('log')

pdat,pobs,pobs_nr,pobs_r = get_data.get_data('P')
plotit(pdat,pobs,ax,
    'DarkBlue',#'#0099CC',
    'o','Praesepe','ADAMSPT','NUM_SPECTRA',0,'lha_noerrs')
hdat,hobs,hobs_nr,hobs_r = get_data.get_data('H')
plotit(hdat,hobs,ax,
    'OrangeRed',#'#FFAA00',
    'D','Hyades','MDM_SPEC_ADAMSPT','MDM_SPECMATCH',0,'lha_noerrs')
ax.legend(numpoints=1,handletextpad=0.2,
    handlelength=1,borderaxespad=0.2,loc=2)
ax.set_ylim(5e-7,1e-3)
"""
