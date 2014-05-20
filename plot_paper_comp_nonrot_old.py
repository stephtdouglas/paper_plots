import numpy as np
import get_data
import matplotlib.pyplot as plt

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



split_bins = np.array([0.,0.3,0.35,0.4,0.65,1.1])
split_bins = np.arange(0.1,1.15,0.1)
split_colp = pdat.field('KH_MASS')
split_colh = hdat.field('KH_MASS')

bin_centers = np.arange(0.15,1.2,0.1)
avg_lha_bin_rot = np.ones(len(bin_centers))*-99.
unc_lha_bin_rot = np.ones(len(bin_centers))*-99.
avg_eqw_bin_rot = np.ones(len(bin_centers))*-99.
unc_eqw_bin_rot = np.ones(len(bin_centers))*-99.
avg_lha_bin_nrot = np.ones(len(bin_centers))*-99.
unc_lha_bin_nrot = np.ones(len(bin_centers))*-99.
avg_eqw_bin_nrot = np.ones(len(bin_centers))*-99.
unc_eqw_bin_nrot = np.ones(len(bin_centers))*-99.
avg_lha_rot = np.ones(len(bin_centers))*-99.
unc_lha_rot = np.ones(len(bin_centers))*-99.
avg_eqw_rot = np.ones(len(bin_centers))*-99.
unc_eqw_rot = np.ones(len(bin_centers))*-99.
avg_lha_nrot = np.ones(len(bin_centers))*-99.
unc_lha_nrot = np.ones(len(bin_centers))*-99.
avg_eqw_nrot = np.ones(len(bin_centers))*-99.
unc_eqw_nrot = np.ones(len(bin_centers))*-99.

"""
for i in range(len(split_bins)-1):
    py = pll
    puy = pull
    hy = hll
    huy = hull
    psplit_stars = np.where((split_colp>(split_bins[i])) & 
        (split_colp<(split_bins[i+1])) & (py>0))[0]
    hsplit_stars = np.where((split_colh>(split_bins[i])) & 
        (split_colh<(split_bins[i+1])) & (hy>0))[0]
    psplit_obsr = np.intersect1d(pobsr,psplit_stars)
    psplit_obsbr = np.intersect1d(pobsbr,psplit_stars)
    psplit_obsnr = np.intersect1d(pobsnr,psplit_stars)
    psplit_obsbnr = np.intersect1d(pobsbnr,psplit_stars)
    hsplit_obsr = np.intersect1d(hobsr,hsplit_stars)
    hsplit_obsbr = np.intersect1d(hobsbr,hsplit_stars)
    hsplit_obsnr = np.intersect1d(hobsnr,hsplit_stars)
    hsplit_obsbnr = np.intersect1d(hobsbnr,hsplit_stars)

    avg_lha_bin_rot[i] = np.average(np.append(py[psplit_obsbr],
         hy[hsplit_obsbr]))
#    print avg_lha_bin_rot[i]
#    print np.append(py[psplit_obsbr],
#         hy[hsplit_obsbr])
    unc_lha_bin_rot[i] = np.sqrt(np.sum(np.append(py[psplit_obsbr]**2,
         hy[hsplit_obsbr]**2)))/(len(psplit_obsbr)+len(hsplit_obsbr))

    avg_lha_bin_nrot[i] = np.average(np.append(py[psplit_obsbnr],
         hy[hsplit_obsbnr]))
    unc_lha_bin_nrot[i] = np.sqrt(np.sum(np.append(py[psplit_obsbnr]**2,
         hy[hsplit_obsbnr]**2)))/(len(psplit_obsbnr)+len(hsplit_obsbnr))

    avg_lha_rot[i] = np.average(np.append(py[psplit_obsr],hy[hsplit_obsr]))
    unc_lha_rot[i] = np.sqrt(np.sum(np.append(py[psplit_obsr]**2,
         hy[hsplit_obsr]**2)))/(len(psplit_obsr)+len(hsplit_obsr))

    avg_lha_nrot[i] = np.average(np.append(py[psplit_obsnr],
         hy[hsplit_obsnr]))
    unc_lha_nrot[i] = np.sqrt(np.sum(np.append(py[psplit_obsnr]**2,
         hy[hsplit_obsnr]**2)))/(len(psplit_obsnr)+len(hsplit_obsnr))
"""
def add_it(avg1,unc1,avg2,unc2):
    avg_out = copy(avg1)
    unc_out = copy(unc1)
    good = np.where((isnan(avg1)==False) & (isnan(avg2)==False)
        & (avg1>-98) & (avg2>-98))[0]
    avg_out[good] = (avg1[good]+avg2[good])/2.
    unc_out[good] = np.sqrt(unc1[good]**2+unc2[good]**2)/2.
    good2 = np.where(isnan(avg1) & (isnan(avg2)==False)  & (avg2>-98))[0]
    avg_out[good2] = avg2[good]
    unc_out[good2] = unc2[good]
    return avg_out, unc_out
"""
plt.figure(figsize=(12,5))
ax = plt.subplot(121)
ax.errorbar(bin_centers,avg_lha_rot,unc_lha_rot,fmt='o',color='r',mec='r',ms=9,
    label='Rotators')
ax.errorbar(bin_centers,avg_lha_nrot,unc_lha_nrot,fmt='^',color='k',ms=9,
    label='Non-rotators')
avg_lha_rot,unc_lha_rot = add_it(avg_lha_rot,unc_lha_rot,
    avg_lha_bin_rot,unc_lha_bin_rot)
avg_lha_nrot,unc_lha_nrot = add_it(avg_lha_nrot,unc_lha_nrot,
    avg_lha_bin_nrot,unc_lha_bin_nrot)
ax.errorbar(bin_centers[4:]+0.01,avg_lha_rot[4:],unc_lha_rot[4:],fmt='o',
    color='r',mec='r',mfc='None',ms=9,label='Rotators (w/ binaries)')
ax.errorbar(bin_centers[4:]+0.01,avg_lha_nrot[4:],unc_lha_nrot[4:],fmt='^',
    color='k',mfc='None',ms=9,label='Non-rotators (w/ binaries)')
ax.set_xlim(0.7,0.1)
ax.set_ylim(1e-6,1e-3)
ax.set_yscale('log')
ax.set_xlabel(r'$M (M_{Sun})$',fontsize='x-large')
ax.set_ylabel(r'$L_{H\alpha}/L_{bol}$',fontsize='x-large')
ax.tick_params(labelsize='large')
ax.legend(loc='best',numpoints=1,handletextpad=0.2)
"""


for i in range(len(split_bins)-1):
    py = peqw
    puy = pueqw
    hy = heqw
    huy = hueqw
    psplit_stars = np.where((split_colp>(split_bins[i])) & 
        (split_colp<(split_bins[i+1])) & (py>-20))[0]
    hsplit_stars = np.where((split_colh>(split_bins[i])) & 
        (split_colh<(split_bins[i+1])) & (hy>-20))[0]
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


#ax = plt.subplot(122)
plt.figure()
ax = plt.subplot(111)
ax.errorbar(bin_centers,avg_eqw_rot,unc_eqw_rot,fmt='o',color='Lime',mec='Lime',ms=9,
    label='Rotators')
ax.errorbar(bin_centers,avg_eqw_nrot,unc_eqw_nrot,fmt='^',color='k',ms=9,
    label='Non-rotators')
avg_eqw_rot,unc_eqw_rot = add_it(avg_eqw_rot,unc_eqw_rot,
    avg_eqw_bin_rot,unc_eqw_bin_rot)
avg_eqw_nrot,unc_eqw_nrot = add_it(avg_eqw_nrot,unc_eqw_nrot,
    avg_eqw_bin_nrot,unc_eqw_bin_nrot)
ax.errorbar(bin_centers[4:]+0.01,avg_eqw_rot[4:],unc_eqw_rot[4:],fmt='o',
    color='Lime',mec='Lime',mfc='None',ms=9,label='Rotators (w/ binaries)')
ax.errorbar(bin_centers[4:]+0.01,avg_eqw_nrot[4:],unc_eqw_nrot[4:],fmt='^',
    color='k',mfc='None',ms=9,label='Non-rotators (w/ binaries)')
ax.set_xlim(1.1,0.1)
ax.set_ylim(3,-7)
ax.set_xlabel(r'$M (M_{\odot})$',fontsize='x-large')
ax.set_ylabel(r'$H\alpha\ \ EqW$',fontsize='x-large')
ax.tick_params(labelsize='large')
ax.legend(loc=2,numpoints=1,handletextpad=0.2)

#plt.savefig('paper_comp_nonrot.png')
#plt.savefig('paper_comp_nonrot.ps',orientation='landscape')

plt.savefig('paper_eqws_nonrot.png')
plt.savefig('paper_eqws_nonrot.ps')
