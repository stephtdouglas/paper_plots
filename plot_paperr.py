
import pyfits
import asciitable as at
import numpy as np
import convertmass as conv
import get_data
import matplotlib.pyplot as plt

kh = at.read('/home/stephanie/ptf/models/kraushillenbrand5.dat')
ksub = np.array([-23,-19,-15,-12,-11,-9,-7,-6])#arange(-15,-1,1)
kh_rpmK = (kh['Mr'] - 0.035*(kh['Mr']-kh['Mi']) - 0.007 - kh['MK'])[ksub]
kh_spt = kh['SpT'][ksub]
klen = len(kh_spt)

def plot_stars(dat,ax):
    """
    """

    ax.minorticks_on()
    ax.tick_params(which='both',width=1,labelsize='large')


    rmag = dat.field('RPRIME')
    rmK = dat.field('RPRIME_K')
    sourceflag = dat.field('RMAG_FLAG')
    ax.plot(rmK[sourceflag=='T'],rmag[sourceflag=='T'],'*',mfc='Magenta',
        mec='None',label='2MASS USNO/Tycho B,V -> r\'')
    ax.plot(rmK[sourceflag=='C'],rmag[sourceflag=='C'],'o',mfc='#993399',
        mec='None',label='CMC14 r\'')
    ax.plot(rmK[sourceflag=='U'],rmag[sourceflag=='U'],'d',mfc='Orange',
        mec='None',label='UCAC4 r\'')
    ax.plot(rmK[sourceflag=='S'],rmag[sourceflag=='S'],'s',mfc='DarkBlue',
        mec='None',label='SDSS r,i -> r\'')


    ax.set_xlim(0,7)
    ax.set_xlabel('(r\'-K)',fontsize='x-large')
    ax.set_ylabel('r\'',fontsize='x-large')


def sidebyside():
    """
    """
    hdat, hobs, hobs_nr, hobs_r = get_data.get_data('H')
    pdat, pobs, pobs_nr, pobs_r = get_data.get_data('P')

#    plt.figure(figsize=(18,7))
    plt.figure(figsize=(13,5.5))
    axh = plt.subplot(122)
    axp = plt.subplot(121)

    plot_stars(hdat,axh)
    axh.set_ylim(18,3)
    texty = 2.75
    for i in range(klen):
        axh.text(kh_rpmK[i],texty,kh_spt[i],fontsize='large')
    axh.tick_params(which='both',top=False)
    axh.text(5.75,3.75,'Hyades',fontsize='large')

    plot_stars(pdat,axp)
    axp.set_ylim(21,6)
    axp.legend(numpoints=1,prop={'size':12},markerscale=1.5,
         handletextpad=0.3,handlelength=1,borderaxespad=0.2)
    texty = 5.75
    for i in range(klen):
       axp.text(kh_rpmK[i],texty,kh_spt[i],fontsize='large')
    axp.tick_params(which='both',top=False)
    axp.text(5.4,10.4,'Praesepe',fontsize='large')

    plt.savefig('paper_rsource.eps',orientation='landscape',bbox_inches='tight')
    plt.savefig('paper_rsource.png')

sidebyside()
