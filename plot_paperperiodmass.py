import get_data,convertmass
import matplotlib.pyplot as plt
import numpy as np


pdat,pobs,pobsnr,pobsr = get_data.get_data('P')
hdat,hobs,hobsnr,hobsr = get_data.get_data('H')

pbinary = ((pdat.field('RPRIME_K')<=4.) & (pdat.field('BINARY')>0))
hbinary = ((hdat.field('RPRIME_K')<=4.) & (hdat.field('BINARY')>0))

pmass = pdat.field('KH_MASS')
hmass = hdat.field('KH_MASS')

kun = hdat.field('KUNDERT_PROT')
delh = hdat.field('DELORME_LITP')
delp = pdat.field('SWASP_PERIOD')
sch = pdat.field('SCHOLZ_PERIOD')
ptf = pdat.field('PTF_PERIOD')
pperiods = pdat.field('PERIOD')
hperiods = hdat.field('PERIOD')
pflag = pdat.field('PERIOD_FLAG')
hflag = hdat.field('PERIOD_FLAG')
ptf_flag = pdat.field('PTF_FLAG')

bad = np.where((kun>0) & (delh>0) & (abs(kun-delh)>0.1))[0]
print bad
delh[bad] = -99.

plt.figure(figsize=(9,8))
ax = plt.subplot(111)
"""
ax.plot(hmass[hflag=='D'],delh[hflag=='D'],'D',mec='#00CCFF',mfc='None')
ax.plot(hmass,kun,'D',mec='Green',mfc='None')
ax.plot(pmass,delp,'o',mec='#00CCFF',mfc='None')
ax.plot(pmass,sch,'o',mec='SlateBlue',mfc='None')
ax.plot(pmass[ptf_flag!='P2m'],ptf[ptf_flag!='P2m'],'o',mec='Red',mfc='None')

ax.plot(hmass[hbinary==False],kun[hbinary==False],'D',mec='Green',
    mfc='DarkGreen')
ax.plot(hmass[(hflag=='D') & (hbinary==False)],delh[(hflag=='D') & (hbinary==False)],'D',mec='#00CCFF',
    mfc='#00CCFF')
ax.plot(pmass[pbinary==False],delp[pbinary==False],'o',mec='#00CCFF',
    mfc='#00CCFF')
ax.plot(pmass[pbinary==False],sch[pbinary==False],'o',mec='SlateBlue',
    mfc='DarkSlateBlue')
ax.plot(pmass[(ptf_flag!='P2m') & (pbinary==False)],
    ptf[(ptf_flag!='P2m') &(pbinary==False)],'o',mec='Red',
    mfc='Red')
"""
ax.plot(hmass[hflag=='D'],hperiods[hflag=='D'],
    'D',mec='#00CCFF',mfc='None')
ax.plot(hmass[hflag=='K'],hperiods[hflag=='K'],
    'D',mec='Green',mfc='None')
ax.plot(hmass[hflag=='R'],hperiods[hflag=='R'],
    'D',mec='Goldenrod',mfc='None')
ax.plot(pmass[pflag=='D'],pperiods[pflag=='D'],
    'o',mec='#00CCFF',mfc='None')
ax.plot(pmass[pflag=='P'],pperiods[pflag=='P'],
    'o',mec='r',mfc='None')
ax.plot(pmass[pflag=='S'],pperiods[pflag=='S'],
    'o',mec='SlateBlue',mfc='None')


ax.plot(hmass[(hflag=='D') & (hbinary==False)],
    hperiods[(hflag=='D') & (hbinary==False)],
    'D',mec='#00CCFF',mfc='#00CCFF')
ax.plot(hmass[(hflag=='K') & (hbinary==False)],
    hperiods[(hflag=='K') & (hbinary==False)],
    'D',mec='Green',mfc='Green')
ax.plot(hmass[(hflag=='R') & (hbinary==False)],
    hperiods[(hflag=='R') & (hbinary==False)],
    'D',mec='Goldenrod',mfc='Goldenrod')
ax.plot(pmass[(pflag=='D') & (pbinary==False)],
    pperiods[(pflag=='D') & (pbinary==False)],
    'o',mec='#00CCFF',mfc='#00CCFF')
ax.plot(pmass[(pflag=='P') & (pbinary==False)],
    pperiods[(pflag=='P') & (pbinary==False)],
    'o',mec='r',mfc='r')
ax.plot(pmass[(pflag=='S') & (pbinary==False)],
    pperiods[(pflag=='S') & (pbinary==False)],
    'o',mec='SlateBlue',mfc='SlateBlue')


ax.set_xlim(1.4,0.1)
ax.set_ylim(0.1,50)
ax.set_yscale('log')
ax.set_xlabel(r'Mass ($M_{\odot}$)',fontsize='x-large')
ax.set_ylabel('Period (d)',fontsize='x-large')
ax.tick_params(labelsize='large')
ax.plot(1.33,40,'ko')
ax.plot(1.29,40,'ko',mfc='None')
ax.plot(1.33,28,'kD')
ax.plot(1.29,28,'kD',mfc='None')
ax.text(1.26,37,'Praesepe (Potential Binary)',color='k',fontsize='large')
ax.text(1.26,26,'Hyades (Potential Binary)',color='k',fontsize='large')

texty = np.logspace(np.log10(0.4),np.log10(1.4),num=6)
ax.text(1.35,texty[5],'Periods measured by',color='k',fontsize='large')
ax.text(1.35,texty[4],'Agueros+ 2011',color='r',fontsize='large')
ax.plot((1.294,1.285),(1.26,1.26),'rs',ms=1,mec='r')
ax.text(1.35,texty[3],'Kundert+ in prep',color='Green',fontsize='large')
ax.text(1.35,texty[2],'Delorme+ 2011',color='#0099CC',fontsize='large')
ax.text(1.35,texty[1],'Scholz+ 2007,2011',color='SlateBlue',fontsize='large')
ax.text(1.35,texty[0],'Radick+ 1987,1995',color='Goldenrod',fontsize='large')



def calc_mass(mK,dist,dist_err,mass_err_frac=0.05):
    # dmod = m-M = 5*log10(d) - 5
    # M = m-dmod = m - 5*log10(d) + 5
    # dM/dd = -5/(d*ln(10))
    dmod = 5.0*np.log10(dist) - 5
    MK = mK - dmod

    if type(mass_err_frac)==np.ndarray:
        mass_err_sq = mass_err_frac**2 #not actually a fraction, but real errors
    else:
        mass_err_sq = (mass_err_frac*mK)**2

    dist_err_sq = (dist_err * 5.0 / (dist * np.log(10.0)))**2

    sigma_MK = np.sqrt(mass_err_sq + dist_err_sq)
    print sigma_MK

    MK_low =  MK + sigma_MK #dimmer/lower mass
    MK_high = MK - sigma_MK #brighter/higher mass
    print 'low',MK_low
    print 'high',MK_high

    mass = convertmass.kraus(MK,'K','None')
    mass_low = convertmass.kraus(MK_low,'K','None')
    mass_high = convertmass.kraus(MK_high,'K','None')

    bin_width = 0.1
    mass_bin = np.arange(0.2,1.3,bin_width)
    mass_errs = np.zeros(2*len(mass_bin)).reshape((2,-1))
    print mass_bin

    for i in range(len(mass_bin)):
        m = mass_bin[i]
        m_loc = np.where(abs(mass-m)<(bin_width/2.))[0]
        if len(m_loc)!=1:
            #print m, len(m_loc), 'NOT ONE'
            mass_errs[1][i] = np.average(abs(mass_low-mass)[m_loc])
            mass_errs[0][i] = np.average(abs(mass_high-mass)[m_loc])
        else:
            #print m, len(m_loc)
            #print abs(mass_low-mass)[m_loc], abs(mass_high-mass)[m_loc]
            mass_errs[1][i] = abs(mass_low-mass)[m_loc]
            mass_errs[0][i] = abs(mass_high-mass)[m_loc]


    return mass_bin,mass_errs

def plot_mass(ax,mK,dist,dist_err,mass_err_frac,yval):

    mass_bin,mass_errs = calc_mass(mK,dist,dist_err,mass_err_frac)
    print 'mass errors',mass_errs 

    ax.errorbar(mass_bin,np.ones(len(mass_bin))*yval,xerr=mass_errs,
         lw=0,ecolor='k',color='k',elinewidth=1.5,fmt='o')


dist = 181.5 #pc
dist_err = 6 #pc


mK1 = pdat.field('TWOMASS_K')[pobsr]
mK_loc = mK1.argsort()
mK = mK1[mK_loc]
mK_err1 = pdat.field('TWOMASS_KERR')[pobsr]
mK_err = mK_err1[mK_loc]

plot_mass(ax,mK,dist,dist_err,mK_err,0.12)
ax.text(1.35,0.14,'Typical Mass Uncertainty',color='k')


plt.savefig('paper_periodmass.png')
plt.savefig('paper_periodmass.ps')
