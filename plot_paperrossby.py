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

pgood = ((pbin==False) & (pll>0) 
        & ((ppmem>=pmem_threshold) | (ppmem<0))
        & (pperiods>0))
hgood = ((hbin==False) & (hll>0)  
        & ((hpmem>=pmem_threshold) | (hpmem<0))
        & (hperiods>0))

pgood2 = ((pbin==False) & (pdat.field('NUM_SPECTRA')>0) 
    & ((ppmem>=pmem_threshold) | (ppmem<0)))
hgood2 = ((hbin==False) & (hdat.field('MDM_SPECMATCH')>0))

pgood3 = (pgood2 & (p_ulim_lha<0) & (pll<0))
hgood3 = (hgood2 & (h_ulim_lha<0) & (hll<0))

plow_mass = ((pmass>0) & (pmass<=0.295))
hlow_mass = ((hmass>0) & (hmass<=0.295))

phigh_mass = ((pmass<=1.3) & (pmass>0.295))
hhigh_mass = ((hmass<=1.3) & (hmass>0.295))

pfast = ((p_ppmax<0.75) & (p_ppmax>0))
hfast = ((h_ppmax<0.75) & (h_ppmax>0))

pslow = (p_ppmax>=0.75)
hslow = (h_ppmax>=0.75)

plt.figure(figsize=(9,8))

alph=0.75

ax = subplot(111)
#Plot statistically inactive stars as open symbols
#high mass, slow
ax.plot(pmass[pgood3 & phigh_mass & pslow],
    pperiods[pgood3 & phigh_mass & pslow],'*',
    mec='k',mfc='None',ms=10)
print pmass[pgood3 & phigh_mass & pslow]
ax.plot(hmass[hgood3 & hhigh_mass & hslow],
    hperiods[hgood3 & hhigh_mass & hslow],'*',
    mec='k',mfc='None',ms=10)
print hmass[hgood3 & hhigh_mass & hslow]
#high mass, fast
ax.plot(pmass[pgood3 & phigh_mass & pfast],
    pperiods[pgood3 & phigh_mass & pfast],'*',
    c='k',mfc='None',ms=10,label=r'No $H\alpha$ emission')
print pmass[pgood3 & phigh_mass & pfast]
ax.plot(hmass[hgood3 & hhigh_mass & hfast],
    hperiods[hgood3 & hhigh_mass & hfast],'*',
    c='k',mfc='None',ms=10)
print hmass[hgood3 & hhigh_mass & hfast]

pmarker_sizes = np.floor(np.round(35*(15+2*np.log10(pll))))
hmarker_sizes = np.floor(np.round(35*(15+2*np.log10(hll))))
pmarker_sizes[np.isnan(pmarker_sizes)] = 0
hmarker_sizes[np.isnan(hmarker_sizes)] = 0

print np.where(np.isnan(pmass))
print np.where(np.isnan(hmass))
print np.where(np.isnan(pperiods))
print np.where(np.isnan(hperiods))
print np.where(np.isnan(pll[pll>0]))
print np.where(np.isnan(hll[hll>0]))
print np.where(np.isnan(p_ulim_lha))
print np.where(np.isnan(h_ulim_lha))
print np.where(np.isnan(pmarker_sizes))
print np.where(np.isnan(hmarker_sizes))

#Plot stars with Lha/Lbol>0
#high mass, slow
ax.scatter(pmass[pgood & phigh_mass & pslow],
    pperiods[pgood & phigh_mass & pslow],c='k',marker='*',edgecolors='none',
    s=pmarker_sizes[pgood & phigh_mass & pslow],
    label=r'$H\alpha$ emission')
#print pmarker_sizes[pgood & phigh_mass & pslow]
ax.scatter(hmass[hgood & hhigh_mass & hslow],
    hperiods[hgood & hhigh_mass & hslow],c='k',marker='*',edgecolors='none',
    s=hmarker_sizes[hgood & hhigh_mass & hslow])
#print hmarker_sizes[hgood & hhigh_mass & hslow]
#high mass, fast
ax.scatter(pmass[pgood & phigh_mass & pfast],
    pperiods[pgood & phigh_mass & pfast],c='k',marker='*',edgecolors='none',
    s=pmarker_sizes[pgood & phigh_mass & pfast],
    label=r'$H\alpha$ emission')
#print pmarker_sizes[pgood & phigh_mass & pfast]
ax.scatter(hmass[hgood & hhigh_mass & hfast],
    hperiods[hgood & hhigh_mass & hfast],c='k',marker='*',edgecolors='none',
    s=hmarker_sizes[hgood & hhigh_mass & hfast])
#print hmarker_sizes[hgood & hhigh_mass & hfast]
#low mass
###### PROBLEM LINE
ax.scatter(pmass[pgood & plow_mass],
    pperiods[pgood & plow_mass],c='DarkGrey',marker='*',edgecolors='none',
    s=pmarker_sizes[pgood & plow_mass],
    label=r'$H\alpha$ emission')
#print pmarker_sizes[pgood & plow_mass]
ax.scatter(hmass[hgood & hlow_mass],
    hperiods[hgood & hlow_mass],c='DarkGrey',marker='*',edgecolors='none',
    s=hmarker_sizes[hgood & hlow_mass])
#print hmarker_sizes[hgood & hlow_mass]

#Plot stars with upper limits
#high mass, slow
ax.plot(pmass[pgood2 & (p_ulim_lha>0) & phigh_mass & pslow],
    pperiods[pgood2 & (p_ulim_lha>0) & phigh_mass & pslow],'*',
    mfc='k',mec='none',ms=11)
print pmass[pgood2 & (p_ulim_lha>0) & phigh_mass & pslow]
ax.plot(hmass[hgood2 & (h_ulim_lha>0) & hhigh_mass & hslow],
    hperiods[hgood2 & (h_ulim_lha>0) & hhigh_mass & hslow],'*',
    mfc='k',mec='none',ms=11)
print hmass[hgood2 & (h_ulim_lha>0) & hhigh_mass & hslow]
#high mass, fast
ax.plot(pmass[pgood2 & (p_ulim_lha>0) & phigh_mass & pfast],
    pperiods[pgood2 & (p_ulim_lha>0) & phigh_mass & pfast],'*',
    mfc='k',mec='none',ms=11)
print pmass[pgood2 & (p_ulim_lha>0) & phigh_mass & pfast]
ax.plot(hmass[hgood2 & (h_ulim_lha>0) & hhigh_mass & hfast],
    hperiods[hgood2 & (h_ulim_lha>0) & hhigh_mass & hfast],'*',
    mfc='k',mec='none',ms=11)
print hmass[hgood2 & (h_ulim_lha>0) & hhigh_mass & hfast]

ax.set_yscale('log')
ax.set_ylim(1e-1,40)
ax.set_ylabel('Period (d)',fontsize='x-large')
ax.set_xlabel(r'Mass (M$_\odot$)',fontsize='x-large')
ax.set_xlim(1.15,0.1)
ax.tick_params(labelsize='x-large')

xlims = ax.get_xlim()
mass = np.arange(xlims[1],xlims[0]+0.05,0.05)
logtau = calc_rossby.tau_wright(mass)
tau = 10**logtau
rossby = np.asarray([0.001,0.01,0.1,1])
for r in rossby:
    periods = r*tau
    ax.plot(mass,periods,'k--')
    ax.text(mass[-2],periods[-2]*1.3,'Ro={}'.format(r),fontsize='large')
#    print periods
ax.text(0.3,0.11,'Ro={}'.format(0.001),fontsize='large')


plt.savefig('paper_rossby.png')
plt.savefig('paper_rossby.ps')
