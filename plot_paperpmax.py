import get_data, calc_pmax3
import matplotlib.pyplot as plt


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

bad = np.where((kun>0) & (delh>0) & (abs(kun-delh)>0.1))[0]
print bad
delh[bad] = -99.

plt.figure()
ax = plt.subplot(111)
ax.plot(hmass,kun,'D',mec='DarkGreen',mfc='None')
ax.plot(hmass,delh,'D',mec='#00CCFF',mfc='None')
ax.plot(pmass,delp,'o',mec='#00CCFF',mfc='None')
ax.plot(pmass,sch,'o',mec='DarkSlateBlue',mfc='None')
ax.plot(pmass,ptf,'o',mec='Red',mfc='None')

ax.plot(hmass[hbinary==False],kun[hbinary==False],'D',mec='DarkGreen',
    mfc='DarkGreen')
ax.plot(hmass[hbinary==False],delh[hbinary==False],'D',mec='#00CCFF',
    mfc='#00CCFF')
ax.plot(pmass[pbinary==False],delp[pbinary==False],'o',mec='#00CCFF',
    mfc='#00CCFF')
ax.plot(pmass[pbinary==False],sch[pbinary==False],'o',mec='DarkSlateBlue',
    mfc='DarkSlateBlue')
ax.plot(pmass[pbinary==False],ptf[pbinary==False],'o',mec='Red',
    mfc='Red')

ax.set_xlim(1.5,0)
ax.set_ylim(0.1,50)
ax.set_yscale('log')
ax.set_xlabel(r'Mass ($M_{\odot}$)',fontsize='large')
ax.set_ylabel('Period (d)',fontsize='large')
ax.tick_params(labelsize='large')
ax.plot(1.45,40,'ko')
ax.plot(1.4,40,'ko',mfc='None')
ax.plot(1.45,28,'kD')
ax.plot(1.4,28,'kD',mfc='None')
ax.text(1.37,37,'Praesepe (Potential Binary)',color='k',fontsize='large')
ax.text(1.37,25,'Hyades (Potential Binary)',color='k',fontsize='large')

ax.text(1.45,0.81,'Kundert+ in prep',color='DarkGreen',fontsize='large')
ax.text(1.45,0.6,'Delorme+ 2011',color='#0099CC',fontsize='large')
ax.text(1.45,0.45,'Scholz+ 2007,2011',color='DarkSlateBlue',fontsize='large')
ax.text(1.45,0.33,'Agueros+ 2011',color='r',fontsize='large')
ax.plot((1.377,1.367),(0.4,0.4),'rs',ms=1,mec='r')

outm1,smin1,smax1 = calc_pmax3.interpall_mass()
x1 = np.arange(smin1+0.00001,smax1-0.00005,0.00001)
print min(x1),max(x1)
plt.plot(x1,10**(outm1(x1)),'-',lw=2,color='Grey')
plt.savefig('paper_pmax.png')
plt.savefig('paper_pmax.ps')
