import get_data, calc_pmax3
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
ax.text(1.26,25,'Hyades (Potential Binary)',color='k',fontsize='large')

ax.text(1.35,0.63,'Periods measured by',color='k',fontsize='large')
ax.text(1.35,0.5,'Agueros+ 2011',color='r',fontsize='large')
ax.plot((1.294,1.285),(0.58,0.58),'rs',ms=1,mec='r')
ax.text(1.35,0.395,'Kundert+ in prep',color='Green',fontsize='large')
ax.text(1.35,0.31,'Delorme+ 2011',color='#0099CC',fontsize='large')
ax.text(1.35,0.248,'Scholz+ 2007,2011',color='SlateBlue',fontsize='large')
ax.text(1.35,0.2,'Radick+ 1987,1995',color='Goldenrod',fontsize='large')

outm1,smin1,smax1 = calc_pmax3.interpall_mass()
x1 = np.arange(0.3,smax1-0.00005,0.00001)
#x1 = np.arange(0.33,smax1-0.00005,0.00001)
print min(x1),max(x1)
plt.plot(x1,10**(outm1(x1)),'-',lw=1.5,color='Grey')
plt.plot(x1,10**(outm1(x1))*0.75,'--',lw=1.5,color='Grey')

plt.savefig('paper_periodmass.png')
plt.savefig('paper_periodmass.ps')
