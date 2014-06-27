import numpy as np
import get_data
import matplotlib.pyplot as plt
import asciitable as at

kh = at.read('/home/stephanie/ptf/models/kraushillenbrand5.dat')
ksub = np.array([-23,-19,-15,-12,-11,-9,-7,-6])#arange(-15,-1,1)
kh_rpmK = (kh['Mr'] - 0.035*(kh['Mr']-kh['Mi']) - 0.007 - kh['MK'])[ksub]
kh_spt = kh['SpT'][ksub]
klen = len(kh_spt)

pdat,pobs,pobsnr,pobsr = get_data.get_data('P')
hdat,hobs,hobsnr,hobsr = get_data.get_data('H')


not_binaries = np.where((pdat.field('BINARY')==0) | ((pdat.field('BINARY')==1)
    & (pdat.field('RPRIME_K')>4)))[0]
yes_binaries = np.where((pdat.field('BINARY')==1)
    & (pdat.field('RPRIME_K')<4))[0]

pobsbr = np.intersect1d(pobsr,yes_binaries)
pobsbnr = np.intersect1d(pobsnr,yes_binaries)
pobsr = np.intersect1d(pobsr,not_binaries)
pobsnr = np.intersect1d(pobsnr,not_binaries)

not_binaries = np.where((hdat.field('BINARY')==0) | ((hdat.field('BINARY')==1)
    & (hdat.field('RPRIME_K')>4)))[0]
yes_binaries = np.where((hdat.field('BINARY')==1)
    & (hdat.field('RPRIME_K')<4))[0]

hobsbnr = np.intersect1d(hobsnr,yes_binaries)
hobsbr = np.intersect1d(hobsr,yes_binaries)
hobsr = np.intersect1d(hobsr,not_binaries)
hobsnr = np.intersect1d(hobsnr,not_binaries)

data = {'mdm':'MDM','mdm2':'MDM2','kafka':'KAFKA','sdss':'SDSS','allen':'ALLEN','mage':'MAGE','hydra':'HYDRA'}
dkeys = data.keys()
pperiods = np.ones(len(pdat))*-99.
pperiods[pdat.field('SCHOLZ_PERIOD')>0] = \
    pdat.field('SCHOLZ_PERIOD')[pdat.field('SCHOLZ_PERIOD')>0]
pperiods[pdat.field('SWASP_PERIOD')>0] = \
    pdat.field('SWASP_PERIOD')[pdat.field('SWASP_PERIOD')>0]
pperiods[pdat.field('PTF_PERIOD')>0] = \
    pdat.field('PTF_PERIOD')[pdat.field('PTF_PERIOD')>0]
pkeys = dkeys
prmark = 'o'
pnrmark = '.'
peqw, pueqw = pdat.field('AVG_EQW'),pdat.field('AVG_EQW_ERR')
pll, pull = pdat.field('AVG_LHA'),pdat.field('AVG_LHA_ERR')

hperiods = np.ones(len(hdat))*-99.
hperiods[hdat.field('KUNDERT_PROT')>0] = \
    hdat.field('KUNDERT_PROT')[hdat.field('KUNDERT_PROT')>0]
hperiods[hdat.field('DELORME_LITP')>0] = \
    hdat.field('DELORME_LITP')[hdat.field('DELORME_LITP')>0]
hkeys = ['mdm']
hrmark = 'D'
hnrmark = 'x'
heqw, hueqw = hdat.field('AVG_EQW'),hdat.field('AVG_EQW_ERR')
hll,  hull = hdat.field('AVG_LHA'),hdat.field('AVG_LHA_ERR')

hpmem = np.ones(len(hdat))*-9999
hpmem[hdat.field('ROESER_DISTANCE')<=9.0] = 99.0
hpmem[(hdat.field('ROESER_DISTANCE')>9.0) & 
     (hdat.field('ROESER_DISTANCE')<=18.0)] = 92.5
hpmem[hdat.field('ROESER_DISTANCE')>18.0] = 70.0
ppmem = pdat.field('ADAMPMEM')
pmem_threshold = 70.0


colorbins = np.logspace(np.log10(0.5),np.log10(8),16)
colorcenter  = (colorbins[1:]+colorbins[:-1])/2.0
num_bins = len(colorcenter)
coloredges = np.zeros(num_bins*2).reshape(2,-1)
coloredges[0] = colorbins[1:]-colorcenter
coloredges[1] = colorcenter-colorbins[:-1]

split_bins = colorbins
split_colp = pdat.field('RPRIME_K')
split_colh = hdat.field('RPRIME_K')

bin_centers = colorcenter
avg_eqw_bin_rot = np.ones(len(bin_centers))*-99.
unc_eqw_bin_rot = np.ones(len(bin_centers))*-99.
avg_eqw_bin_nrot = np.ones(len(bin_centers))*-99.
unc_eqw_bin_nrot = np.ones(len(bin_centers))*-99.
avg_eqw_rot = np.ones(len(bin_centers))*-99.
unc_eqw_rot = np.ones(len(bin_centers))*-99.
avg_eqw_nrot = np.ones(len(bin_centers))*-99.
unc_eqw_nrot = np.ones(len(bin_centers))*-99.

def add_it(avg1,unc1,avg2,unc2):
    avg_out = np.copy(avg1)
    unc_out = np.copy(unc1)
    good = np.where((np.isnan(avg1)==False) & (np.isnan(avg2)==False)
        & (avg1>-98) & (avg2>-98))[0]
    avg_out[good] = (avg1[good]+avg2[good])/2.
    unc_out[good] = np.sqrt(unc1[good]**2+unc2[good]**2)/2.
    good2 = np.where(np.isnan(avg1) & (np.isnan(avg2)==False)  & (avg2>-98))[0]
    avg_out[good2] = avg2[good]
    unc_out[good2] = unc2[good]
    return avg_out, unc_out

for i in range(len(split_bins)-1):
    py = peqw
    puy = pueqw
    hy = heqw
    huy = hueqw
    psplit_stars = np.where((split_colp>(split_bins[i])) & 
        (split_colp<(split_bins[i+1])) & (py>-20)
        & ((ppmem>=pmem_threshold) | (ppmem<0)) )[0]
    hsplit_stars = np.where((split_colh>(split_bins[i])) & 
        (split_colh<(split_bins[i+1])) & (hy>-20)
        & ((hpmem>=pmem_threshold) | (hpmem<0)) )[0]
    psplit_obsr = np.intersect1d(pobsr,psplit_stars)
    psplit_obsbr = np.intersect1d(pobsbr,psplit_stars)
    psplit_obsnr = np.intersect1d(pobsnr,psplit_stars)
    psplit_obsbnr = np.intersect1d(pobsbnr,psplit_stars)
    hsplit_obsr = np.intersect1d(hobsr,hsplit_stars)
    hsplit_obsbr = np.intersect1d(hobsbr,hsplit_stars)
    hsplit_obsnr = np.intersect1d(hobsnr,hsplit_stars)
    hsplit_obsbnr = np.intersect1d(hobsbnr,hsplit_stars)

    avg_eqw_bin_rot[i] = -1*np.average(np.append(py[psplit_obsbr],
         hy[hsplit_obsbr]))
    unc_eqw_bin_rot[i] = np.std(np.append(py[psplit_obsbr],
         hy[hsplit_obsbr]))
#    unc_eqw_bin_rot[i] = np.sqrt(np.sum(np.append(puy[psplit_obsbr]**2,
#         huy[hsplit_obsbr]**2)))/(len(psplit_obsbr)+len(hsplit_obsbr))

    avg_eqw_bin_nrot[i] = -1*np.average(np.append(py[psplit_obsbnr],
         hy[hsplit_obsbnr]))
    unc_eqw_bin_nrot[i] = np.std(np.append(py[psplit_obsbnr],
         hy[hsplit_obsbnr]))
#    unc_eqw_bin_nrot[i] = np.sqrt(np.sum(np.append(puy[psplit_obsbnr]**2,
#         huy[hsplit_obsbnr]**2)))/(len(psplit_obsbnr)+len(hsplit_obsbnr))

    avg_eqw_rot[i] = -1*np.average(np.append(py[psplit_obsr],hy[hsplit_obsr]))
    unc_eqw_rot[i] = np.std(np.append(py[psplit_obsr],hy[hsplit_obsr]))
#    unc_eqw_rot[i] = np.sqrt(np.sum(np.append(puy[psplit_obsr]**2,
#         huy[hsplit_obsr]**2)))/(len(psplit_obsr)+len(hsplit_obsr))

    avg_eqw_nrot[i] = -1*np.average(np.append(py[psplit_obsnr],
         hy[hsplit_obsnr]))
    unc_eqw_nrot[i] = np.std(np.append(py[psplit_obsnr],
         hy[hsplit_obsnr]))
#    unc_eqw_nrot[i] = np.sqrt(np.sum(np.append(puy[psplit_obsnr]**2,
#         huy[hsplit_obsnr]**2)))/(len(psplit_obsnr)+len(hsplit_obsnr))

plt.figure(figsize=(9,8))
ax = plt.subplot(111)
ax.errorbar(bin_centers,avg_eqw_rot,unc_eqw_rot,coloredges,
    fmt='*',color='#FF4D4D',mec='#FF4D4D',ms=9,label='Rotators',capsize=0)
ax.errorbar(bin_centers,avg_eqw_nrot,unc_eqw_nrot,coloredges,
    fmt='^',color='k',ms=9,capsize=0,
    label='Non-periodic variables')
avg_eqw_rot,unc_eqw_rot = add_it(avg_eqw_rot,unc_eqw_rot,
    avg_eqw_bin_rot,unc_eqw_bin_rot)
avg_eqw_nrot,unc_eqw_nrot = add_it(avg_eqw_nrot,unc_eqw_nrot,
    avg_eqw_bin_nrot,unc_eqw_bin_nrot)
ax.errorbar(bin_centers[:-3]+0.03,avg_eqw_rot[:-3],unc_eqw_rot[:-3],fmt='*',
    color='#FF4D4D',mec='#FF4D4D',mfc='None',ms=9,
    label='Rotators (w/ binaries)',capsize=0)
ax.errorbar(bin_centers[:-3]+0.03,avg_eqw_nrot[:-3],unc_eqw_nrot[:-3],fmt='^',
    color='k',mfc='None',ms=9,label='Non-periodic variables (w/ binaries)',capsize=0)
ax.set_xlim(0.25,6.8)
ax.set_ylim(8,-10)
ax.set_xlabel('(r\'-K)',fontsize='large')
ax.set_ylabel(r'$H\alpha\ \ EqW$',fontsize='x-large')
ax.legend(loc=2,numpoints=1,handletextpad=0.2)
ax.tick_params(which='both',labelsize='large',top=False)
texty = -10.25
for i in range(klen):
   ax.text(kh_rpmK[i],texty,kh_spt[i],fontsize='large')


plt.savefig('paper_eqws_nonrot.png')
plt.savefig('paper_eqws_nonrot.ps')
