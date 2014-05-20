import numpy as np
import get_data, calc_rossby, calc_pmax2, fitting
import matplotlib.pyplot as plt
import plot_grid as pg
import asciitable as at


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

pdat,pobs,pobsnr,pobsr = get_data.get_data('P')
hdat,hobs,hobsnr,hobsr = get_data.get_data('H')
plen = len(pdat)
hlen = len(hdat)
pra = pdat.field('RA')
pdec = pdat.field('DEC')
hra = hdat.field('RA')
hdec = hdat.field('DEC')
pbin = (pdat.field('BINARY')>0)
hbin = (hdat.field('BINARY')>0)
pperiods = pdat.field('PERIOD')
hperiods = hdat.field('PERIOD')
pmass = pdat.field('KH_MASS')
hmass = hdat.field('KH_MASS')
p_ppmax = calc_pmax2.p_pmax(pmass,pperiods)
h_ppmax = calc_pmax2.p_pmax(hmass,hperiods)

hpmem = hdat.field('ROESER_PMEM')
ppmem = pdat.field('ADAMPMEM')
pmem_threshold = 70.0


rpmK = np.append(pdat.field('RPRIME_K'),hdat.field('RPRIME_K'))
rpmKerr = np.append(pdat.field('RPRIME_K_ERR'),hdat.field('RPRIME_K_ERR'))
peqw,pueqw = pdat.field('AVG_EQW'),pdat.field('AVG_EQW_ERR')
pll,pull = pdat.field('AVG_LHA'),pdat.field('AVG_LHA_ERR')
heqw,hueqw = hdat.field('AVG_EQW'),hdat.field('AVG_EQW_ERR')
hll,hull = hdat.field('AVG_LHA'),hdat.field('AVG_LHA_ERR')
eqw = np.append(peqw,heqw)
ueqw = np.append(pueqw,hueqw)
binary = np.append(pbin,hbin)
chi = np.append(pdat.field('CHI'),hdat.field('CHI'))

pros = pdat.field('ROSSBY')
hros = hdat.field('ROSSBY')


p_ulim_lha, p_err_ulim_lha, pdnerr = simple(peqw,pueqw,pdat.field('CHI'))
h_ulim_lha, h_err_ulim_lha, hdnerr = simple(heqw,hueqw,hdat.field('CHI'))

pmem_threshold = 70.0

plt.figure(figsize=(13.5,6))
ax = subplot(122)
ax.plot((0.13,0.13),(1e-7,1e-3),'k--')
ax.errorbar(
    10**pros[(p_ppmax<0.75) & (pmass<=1.3) & (pmass>0.3) & (pbin==False) 
        & ((ppmem>=pmem_threshold) | (ppmem<0))],
    pll[(p_ppmax<0.75) & (pmass<=1.3) & (pmass>0.3) & (pbin==False) 
        & ((ppmem>=pmem_threshold) | (ppmem<0))],
    pull[(p_ppmax<0.75) & (pmass<=1.3) & (pmass>0.3) & (pbin==False) 
        & ((ppmem>=pmem_threshold) | (ppmem<0))],
    color='k',fmt='*',capsize=0,ms=11,mec='k')
ax.errorbar(
    10**hros[(h_ppmax<0.75) & (hmass<=1.3) & (hmass>0.3) & (hbin==False)  
        & ((hpmem>=pmem_threshold) | (hpmem<0))],
    hll[(h_ppmax<0.75) & (hmass<=1.3) & (hmass>0.3) & (hbin==False) 
        & ((hpmem>=pmem_threshold) | (hpmem<0))],
    hull[(h_ppmax<0.75) & (hmass<=1.3) & (hmass>0.3) & (hbin==False) 
        & ((hpmem>=pmem_threshold) | (hpmem<0))],
    color='k',fmt='*',capsize=0,ms=11,mec='k')
ax.errorbar(10**pros[(pmass<=0.3) & ((ppmem>=pmem_threshold) | (ppmem<0))],
    pll[(pmass<=0.3) & ((ppmem>=pmem_threshold) | (ppmem<0))],
    pull[(pmass<=0.3) & ((ppmem>=pmem_threshold) | (ppmem<0))],ms=11,
    color='DarkGrey',fmt='*',mfc='DarkGrey',capsize=0,mec='DarkGrey')
ax.errorbar(10**hros[(hmass<=0.3)& ((hpmem>=pmem_threshold) | (hpmem<0))],
    hll[(hmass<=0.3)& ((hpmem>=pmem_threshold) | (hpmem<0))],
    hull[(hmass<=0.3)& ((hpmem>=pmem_threshold) | (hpmem<0))],ms=11,
    color='DarkGrey',fmt='*',mfc='DarkGrey',capsize=0,mec='DarkGrey')
ax.errorbar(
    10**pros[(p_ppmax>0.75) & (pmass<=1.3) & (pmass>0.3) & (pbin==False) 
        & ((peqw-pueqw)>0) & ((ppmem>=pmem_threshold) | (ppmem<0))],
    pll[(p_ppmax>0.75) & (pmass<=1.3) & (pmass>0.3) & (pbin==False) 
        & ((peqw-pueqw)>0) & ((ppmem>=pmem_threshold) | (ppmem<0))],
    pull[(p_ppmax>0.75) & (pmass<=1.3) & (pmass>0.3) & (pbin==False) 
        & ((peqw-pueqw)>0) & ((ppmem>=pmem_threshold) | (ppmem<0))],ms=11,
    color='#FF4D4D',fmt='*',capsize=0,mec='#FF4D4D')
ax.errorbar(
    10**hros[(h_ppmax>0.75) & (hmass<=1.3) & (hmass>0.3) & (hbin==False) 
        & ((heqw-hueqw)>0)& ((hpmem>=pmem_threshold) | (hpmem<0))],
    hll[(h_ppmax>0.75) & (hmass<=1.3) & (hmass>0.3) & (hbin==False) 
        & ((heqw-hueqw)>0)& ((hpmem>=pmem_threshold) | (hpmem<0))],
    hull[(h_ppmax>0.75) & (hmass<=1.3) & (hmass>0.3) & (hbin==False) 
        & ((heqw-hueqw)>0)& ((hpmem>=pmem_threshold) | (hpmem<0))],ms=11,
        color='#FF4D4D',fmt='*',capsize=0,mec='#FF4D4D')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylim(5e-7,6e-4)
ax.set_ylabel(r'$L_{H\alpha}/L_{bol}$',fontsize='xx-large')
ax.set_xlabel('Ro',fontsize='x-large')
ax.set_xlim(1e-3,1)
ax.tick_params(labelsize='x-large')
ax.set_xticklabels((0.001,0.01,0.1,1))
ax.text(0.025,4.25e-4,'Saturated',fontsize='large')
ax.text(0.16,4.25e-4,'Unsaturated',fontsize='large')

#Upper limits turn out to be only for stars that have spun down

pgood = np.where((p_ppmax>0.75) & (pmass<=1.3) & (pmass>0.3) & (pbin==False))[0]
hgood = np.where((h_ppmax>0.75) & (hmass<=1.3) & (hmass>0.3) & (hbin==False))[0]

pup = (np.ones(2*len(pgood))*1e-6).reshape(2,-1)
pup[1] = p_err_ulim_lha[pgood]
pup[0] = pdnerr[pgood]
#print (p_ulim_lha-pdnerr)[pgood]
hup = (np.ones(2*len(hgood))*1e-6).reshape(2,-1)
hup[1] = h_err_ulim_lha[hgood]
hup[0] = hdnerr[hgood] 
#print (h_ulim_lha-hdnerr)[hgood]
#print ppmem[pgood]
#print hpmem[hgood]

ax.errorbar(
    10**pros[pgood],
    p_ulim_lha[pgood],
    pup,ecolor='#FF4D4D',
    color='#FF4D4D',fmt=None,lolims=True,ms=3)
ax.errorbar(
    10**hros[hgood],
    h_ulim_lha[hgood],
    hup,ecolor='#FF4D4D',
    color='#FF4D4D',fmt=None,lolims=True,ms=3)

pgood2 = np.where((p_ppmax<0.75) & (pmass<=1.3) & (pmass>0.3) & (pbin==False))[0]
hgood2 = np.where((h_ppmax<0.75) & (hmass<=1.3) & (hmass>0.3) & (hbin==False))[0]

pup = (np.ones(2*len(pgood2))*1e-6).reshape(2,-1)
pup[1] = p_err_ulim_lha[pgood2]
pup[0] = pdnerr[pgood2]
#print (p_ulim_lha-pdnerr)[pgood2]
hup = (np.ones(2*len(hgood2))*1e-6).reshape(2,-1)
hup[1] = h_err_ulim_lha[hgood2]
hup[0] = hdnerr[hgood2] 
#print (h_ulim_lha-hdnerr)[hgood2]
#print ppmem[pgood2]
#print hpmem[hgood2]

ax.errorbar(
    10**pros[pgood2],
    p_ulim_lha[pgood2],
    pup,ecolor='k',
    color='k',fmt=None,lolims=True,ms=3)
ax.errorbar(
    10**hros[hgood2],
    h_ulim_lha[hgood2],
    hup,ecolor='k',
    color='k',fmt=None,lolims=True,ms=3)






ax = subplot(121)
ax.plot(pmass[(p_ppmax<0.75) & (pmass<=1.3) & (pmass>0.3) & (pbin==False) 
        & ((ppmem>=pmem_threshold) | (ppmem<0))],
    pperiods[(p_ppmax<0.75) & (pmass<=1.3) & (pmass>0.3) & (pbin==False) 
        & ((ppmem>=pmem_threshold) | (ppmem<0))],
    'k*',mfc='None',ms=11,label=r'No $H\alpha$ emission')
ax.plot(hmass[(h_ppmax<0.75) & (hmass<=1.3) & (hmass>0.3) & (hbin==False)  
        & ((hpmem>=pmem_threshold) | (hpmem<0))],
    hperiods[(h_ppmax<0.75) & (hmass<=1.3) & (hmass>0.3) & (hbin==False)  
        & ((hpmem>=pmem_threshold) | (hpmem<0))],
    'k*',mfc='None',ms=11)
ax.plot(pmass[(pmass<=0.3) & (pdat.field('NUM_SPECTRA')>0) 
        & ((ppmem>=pmem_threshold) | (ppmem<0))],
    pperiods[(pmass<=0.3) & (pdat.field('NUM_SPECTRA')>0) 
        & ((ppmem>=pmem_threshold) | (ppmem<0))],
    '*',mfc='None',mec='DarkGrey',ms=11)
ax.plot(hmass[(hmass<=0.3)  
        & ((hpmem>=pmem_threshold) | (hpmem<0))],
    hperiods[(hmass<=0.3)  
        & ((hpmem>=pmem_threshold) | (hpmem<0))],
    '*',mfc='None',mec='DarkGrey',ms=11)
ax.plot(pmass[(p_ppmax>0.75) & (pmass<=1.3) & (pmass>0.3) & (pbin==False) 
        & ((ppmem>=pmem_threshold) | (ppmem<0))],
    pperiods[(p_ppmax>0.75) & (pmass<=1.3) & (pmass>0.3) & (pbin==False) 
        & ((ppmem>=pmem_threshold) | (ppmem<0))],
    '*',mec='#FF4D4D',mfc='None',ms=11)
ax.plot(hmass[(h_ppmax>0.75) & (hmass<=1.3) & (hmass>0.3) & (hbin==False)  
        & ((hpmem>=pmem_threshold) | (hpmem<0))],
    hperiods[(h_ppmax>0.75) & (hmass<=1.3) & (hmass>0.3) & (hbin==False)  
        & ((hpmem>=pmem_threshold) | (hpmem<0))],
    '*',mec='#FF4D4D',mfc='None',ms=11)
ax.plot(pmass[pgood],
    pperiods[pgood],
    '*',mec='#FF4D4D',mfc='None',ms=11)
ax.plot(hmass[hgood],
    hperiods[hgood],
    '*',mec='#FF4D4D',mfc='None',ms=11)
ax.plot(pmass[pgood2],
    pperiods[pgood2],
    '*',mec='k',mfc='None',ms=11)
ax.plot(hmass[hgood2],
    hperiods[hgood2],
    '*',mec='k',mfc='None',ms=11)



print ppmem[(p_ppmax<0.75) & (pmass<=1.3) & (pmass>0.3) & (pbin==False) & (pll<0) 
        & (p_ppmax>0) & ((ppmem>=pmem_threshold) | (ppmem<0))]
ax.plot(pmass[(p_ppmax<0.75) & (pmass<=1.3) & (pmass>0.3) & (pbin==False) & (pll>=0) 
        & ((ppmem>=pmem_threshold) | (ppmem<0))],
    pperiods[(p_ppmax<0.75) & (pmass<=1.3) & (pmass>0.3) & (pbin==False) & (pll>=0) 
        & ((ppmem>=pmem_threshold) | (ppmem<0))],
    'k*',ms=11,label=r'$H\alpha$ emission')
ax.plot(hmass[(h_ppmax<0.75) & (hmass<=1.3) & (hmass>0.3) & (hbin==False) & (hll>=0)  
        & ((hpmem>=pmem_threshold) | (hpmem<0))],
    hperiods[(h_ppmax<0.75) & (hmass<=1.3) & (hmass>0.3) & (hbin==False) & (hll>=0)  
        & ((hpmem>=pmem_threshold) | (hpmem<0))],
    'k*',ms=11)
ax.plot(pmass[(pmass<=0.3) & (pll>=0) & (pdat.field('NUM_SPECTRA')>0) 
        & ((ppmem>=pmem_threshold) | (ppmem<0))],
    pperiods[(pmass<=0.3) & (pll>=0) & (pdat.field('NUM_SPECTRA')>0) 
        & ((ppmem>=pmem_threshold) | (ppmem<0))],
    '*',mfc='DarkGrey',mec='DarkGrey',ms=11)
ax.plot(hmass[(hmass<=0.3) & (hll>=0)  
        & ((hpmem>=pmem_threshold) | (hpmem<0))],
    hperiods[(hmass<=0.3) & (hll>=0)  
        & ((hpmem>=pmem_threshold) | (hpmem<0))],
    '*',mfc='DarkGrey',mec='DarkGrey',ms=11)
ax.plot(pmass[(p_ppmax>0.75) & (pmass<=1.3) & (pmass>0.3) & (pbin==False) & (pll>=0) 
        & ((ppmem>=pmem_threshold) | (ppmem<0))],
    pperiods[(p_ppmax>0.75) & (pmass<=1.3) & (pmass>0.3) & (pbin==False) & (pll>=0) 
        & ((ppmem>=pmem_threshold) | (ppmem<0))],
    '*',color='#FF4D4D',mfc='#FF4D4D',mec='#FF4D4D',ms=11)
ax.plot(hmass[(h_ppmax>0.75) & (hmass<=1.3) & (hmass>0.3) & (hbin==False) & (hll>=0)  
        & ((hpmem>=pmem_threshold) | (hpmem<0))],
    hperiods[(h_ppmax>0.75) & (hmass<=1.3) & (hmass>0.3) & (hbin==False) & (hll>=0)  
        & ((hpmem>=pmem_threshold) | (hpmem<0))],
    '*',color='#FF4D4D',mfc='#FF4D4D',mec='#FF4D4D',ms=11)
ax.plot(pmass[np.intersect1d(pgood,np.where(p_ulim_lha>0)[0])],
    pperiods[np.intersect1d(pgood,np.where(p_ulim_lha>0)[0])],
    '*',mec='#FF4D4D',mfc='#FF4D4D',ms=11)
ax.plot(hmass[np.intersect1d(hgood,np.where(h_ulim_lha>0)[0])],
    hperiods[np.intersect1d(hgood,np.where(h_ulim_lha>0)[0])],
    '*',mec='#FF4D4D',mfc='#FF4D4D',ms=11)
ax.plot(pmass[np.intersect1d(pgood2,np.where(p_ulim_lha>0)[0])],
    pperiods[np.intersect1d(pgood2,np.where(p_ulim_lha>0)[0])],
    'k*',mec='k',ms=11)
ax.plot(hmass[np.intersect1d(hgood2,np.where(h_ulim_lha>0)[0])],
    hperiods[np.intersect1d(hgood2,np.where(h_ulim_lha>0)[0])],
    'k*',mec='k',ms=11)

ax.set_yscale('log')
ax.set_ylim(1e-1,30)
ax.set_ylabel('Period (d)',fontsize='x-large')
ax.set_xlabel(r'Mass (M$_\odot$)',fontsize='x-large')
ax.set_xlim(1.15,0.1)
ax.tick_params(labelsize='x-large')

ax.legend(loc=3,numpoints=1,handletextpad=0.02,frameon=False)

textx = 1.1
ax.text(textx,7e-1,r'Slow Rotators, M>0.3 M$_\odot$',color='#FF4D4D',fontsize='large')
ax.text(textx,5e-1,r'Rapid Rotators, M>0.3 M$_\odot$',color='k',fontsize='large')
ax.text(textx,3.5e-1,r'M<0.3 M$_\odot$',color='DarkGrey',fontsize='large')
ax.text(textx,3.5e-1,r'M<0.3 M$_\odot$',color='DarkGrey',fontsize='large')
ax.plot((1.052,1.03),(0.355,0.355),'-',color='DarkGrey')
ax.plot((1.052,1.03),(0.355,0.355),'-',color='DarkGrey')




plt.savefig('paper_rossby.png')
plt.savefig('paper_rossby.ps',orientation='landscape')
