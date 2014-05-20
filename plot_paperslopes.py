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

ppmem = pdat.field('ADAMPMEM')
hpmem = hdat.field('ROESER_PMEM')
pmem_threshold=70.0


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

max_mass = max(np.append(pmass[(pbin==False) & (peqw-pueqw>0)
        & ((ppmem>=pmem_threshold) | (ppmem<0))],
        hmass[(hbin==False)  & (heqw-hueqw>0)
        & ((hpmem>=pmem_threshold) | (hpmem<0))]))
print max_mass

figure(figsize=(8,12))
ax = subplot(211)
ax.errorbar(
    10**pros[(pmass<=1.3) & (pmass>0.1) & (pbin==False) & (peqw-pueqw>0)
        & ((ppmem>=pmem_threshold) | (ppmem<0))],
    pll[(pmass<=1.3) & (pmass>0.1) & (pbin==False)  & (peqw-pueqw>0)
        & ((ppmem>=pmem_threshold) | (ppmem<0))],
    pull[(pmass<=1.3) & (pmass>0.1) & (pbin==False)  & (peqw-pueqw>0)
        & ((ppmem>=pmem_threshold) | (ppmem<0))],
    color='c',fmt='*',capsize=0,ms=12,mec='c')
ax.errorbar(
    10**hros[(hmass<=1.3) & (hmass>0.1) & (hbin==False)  & (heqw-hueqw>0)
        & ((hpmem>=pmem_threshold) | (hpmem<0))],
    hll[(hmass<=1.3) & (hmass>0.1) & (hbin==False)  & (heqw-hueqw>0)
        & ((hpmem>=pmem_threshold) | (hpmem<0))],
    hull[(hmass<=1.3) & (hmass>0.1) & (hbin==False)  & (heqw-hueqw>0)
        & ((hpmem>=pmem_threshold) | (hpmem<0))],
    color='c',fmt='*',capsize=0,ms=12,mec='c')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylim(7e-7,5e-4)
ax.set_ylabel(r'$L_{H\alpha}/L_{bol}$',fontsize='xx-large')
#ax.set_xlabel('Ro',fontsize='x-large')
ax.set_xlim(1e-3,2)
ax.tick_params(labelsize='x-large')
ax.set_xticklabels((0.001,0.01,0.1,1))

pline = np.where((10**pros>0.001) & (pbin==False) 
        & ((peqw-pueqw)>0) & (10**pros<0.13) 
        & ((ppmem>=pmem_threshold) | (ppmem<0)))[0]
hline = np.where((10**hros>0.001) & (hbin==False) 
        & ((heqw-hueqw)>0) & (10**hros<0.13)  
        & ((hpmem>=pmem_threshold) | (hpmem<0)))[0]
flatlha = np.append(pll[pline],hll[hline])
flatulha = np.append(pull[pline],hull[hline])
sat_level = np.average(flatlha)
plot((0.001,0.13),(sat_level,sat_level),'k-',lw=2)
x = np.arange(0.13,1.2,0.001)
constant = sat_level/(0.13**-1.)
plot(x,constant*(x**-1.),'k-',lw=2,label=r'$\beta=\ -1$')
constant = sat_level/(0.13**-2.1)
plot(x,constant*(x**-2.1),'k--',lw=2,label=r'$\beta=\ -2.1$')
constant = sat_level/(0.13**-2.7)
plot(x,constant*(x**-2.7),'k-.',lw=2,label=r'$\beta=\ -2.7$')
legend(loc=3,title=r'$L_{H\alpha}/L_{bol}\ \propto\ Ro^{\beta}$')
#plt.savefig('uplims_split_{}.png'.format(date))

#xrays
wprae = at.read('../xray/wright_praesepe_catalog.csv')
wra = wprae['ra']
wdec = wprae['dec']
wlen = len(wprae)
pmatch_loch = np.ones(wlen,int)*-99
pmatch_locl = np.ones(wlen,int)*-99
pmatch_h = np.ones(wlen,int)*-99
pmatch_l = np.ones(wlen,int)*-99
wpros = 10**(calc_rossby.rossby(wprae['Prot'],wprae['mass']))
match_tol = 5./3500.
match_count = 0
for p in range(plen):
    for w in range(wlen):
        if ((abs(wra[w]-pra[p])<match_tol) and (abs(wdec[w]-pdec[p])<match_tol)
            and ((ppmem[p]>=pmem_threshold) | (ppmem[p]<0)) 
            & (pbin[p]==False)):
            match_count += 1
            if pmass[p]>max_mass:
                pmatch_loch[w] = p
                pmatch_h[w] = w
            else:
                pmatch_locl[w] = p
                pmatch_l[w] = w
print 'Praesepe:',match_count
pmatch_loch = np.delete(pmatch_loch,np.where(pmatch_h<0))
pmatch_locl = np.delete(pmatch_locl,np.where(pmatch_l<0))
pmatch_h = np.delete(pmatch_h,np.where(pmatch_h<0))
pmatch_l = np.delete(pmatch_l,np.where(pmatch_l<0))

wh = at.read('../xray/wright_hyades_catalog.csv')
wra = wh['ra']
wdec = wh['dec']
wlen = len(wh)
hmatch_loch = np.ones(wlen,int)*-99
hmatch_locl = np.ones(wlen,int)*-99
hmatch_h = np.ones(wlen,int)*-99
hmatch_l = np.ones(wlen,int)*-99
whros = 10**(calc_rossby.rossby(wh['Prot'],wh['mass']))

match_tol = 5./3500.
match_count = 0
rot_match = 0
for h in range(hlen):
    for w in range(wlen):
        if ((abs(wra[w]-hra[h])<match_tol) and (abs(wdec[w]-hdec[h])<match_tol)
            and ((hpmem[h]>=pmem_threshold) | (hpmem[h]<0)) 
            & (hbin[h]==False)):
            match_count += 1
            if hmass[h]>max_mass:
                hmatch_loch[w] = h
                hmatch_h[w] = w
            else:
                hmatch_locl[w] = h
                hmatch_l[w] = w
print 'Hyades:',match_count
hmatch_loch = np.delete(hmatch_loch,np.where(hmatch_h<0))
hmatch_locl = np.delete(hmatch_locl,np.where(hmatch_l<0))
hmatch_h = np.delete(hmatch_h,np.where(hmatch_h<0))
hmatch_l = np.delete(hmatch_l,np.where(hmatch_l<0))

wpLx = 10**wprae['LxLbol']
whLx = 10**wh['LxLbol']

ax = subplot(212)
ax.plot(10**pros[pmatch_loch],wpLx[pmatch_h],'*',ms=12,mec='c',
    mfc='None')
ax.plot(10**hros[hmatch_loch],whLx[hmatch_h],'*',ms=12,mec='c',
    mfc='None')
ax.plot(10**pros[pmatch_locl],wpLx[pmatch_l],'c*',ms=12,mec='c')
ax.plot(10**hros[hmatch_locl],whLx[hmatch_l],'c*',ms=12,mec='c')

sat_level = 10**-3.13
constant = sat_level/(0.13**-2)
ax.plot((0.001,0.13),(sat_level,sat_level),'k-',lw=2)
x = np.arange(0.13,1.2,0.001)
constant = sat_level/(0.13**-1.)
ax.plot(x,constant*(x**-1.),'k-',lw=2,label=r'$\beta=\ -1$')
constant = sat_level/(0.13**-2.1)
ax.plot(x,constant*(x**-2.1),'k--',lw=2,label=r'$\beta=\ -2.1$')
constant = sat_level/(0.13**-2.7)
ax.plot(x,constant*(x**-2.7),'k-.',lw=2,label=r'$\beta=\ -2.7$')
#ax.legend(loc=3)
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_xlim(1e-3,2)
ax.set_ylim(6e-6,5e-3)
ax.set_ylabel(r'$L_X/L_{bol}$',fontsize='xx-large')
ax.set_xlabel('Ro',fontsize='x-large')
ax.tick_params(labelsize='x-large')
ax.set_xticklabels((0.001,0.01,0.1,1))
legend(loc=3,title=r'$L_{X}/L_{bol}\ \propto\ Ro^{\beta}$')

plt.savefig('paperslopes.png')
plt.savefig('paperslopes.ps')


