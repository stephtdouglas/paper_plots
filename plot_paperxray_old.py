import get_data
import asciitable as at
import plot_grid as pg
import calc_pmax, calc_chi, pyfits


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






dkeys = ['mdm','mdm2','allen','kafka','hydra','mage','sdss']

pdat,pobs,pobsnr,pobsr = get_data.get_data('P')
hdat,hobs,hobsnr,hobsr = get_data.get_data('H')
plen=len(pdat)
pra = pdat.field('RA')
pdec = pdat.field('DEC')
peqw,pueqw = pdat.field('AVG_EQW'),pdat.field('AVG_EQW_ERR')
pll,pull = pdat.field('AVG_LHA'),pdat.field('AVG_LHA_ERR')
pperiods = pdat.field('PERIOD')
p_ppmax = calc_pmax.p_pmax(pdat.field('KH_MASS'),pperiods)

pbin = (pdat.field('BINARY')>0)
ppmem = pdat.field('ADAMPMEM')
pmem_threshold=70.0


wprae = at.read('../xray/wright_praesepe_catalog.csv')
wra = wprae['ra']
wdec = wprae['dec']
wlen = len(wprae)
pmatch_loc = np.ones(wlen,int)*-99
pmatch_mloc = np.ones(wlen,int)*-99
pmatch_fgkloc = np.ones(wlen,int)*-99
pmatches_m = np.zeros(wlen,int)
pmatches_fgk = np.zeros(wlen,int)

match_tol = 5./3500.
match_count = 0
for p in range(plen):
    for w in range(wlen):
        if ((abs(wra[w]-pra[p])<match_tol) and (abs(wdec[w]-pdec[p])<match_tol)
            and ((ppmem[p]>=pmem_threshold) or (ppmem[p]<0)) 
            and (pbin[p]==False)
            ):
            #print 'found!',pdat.field('LITNAME')[p],pdat.field('KH_MASS')[p],
            #print pdat.field('NUM_SPECTRA')[p]
            match_count += 1
            pmatch_loc[w] = p
            if pdat.field('ADAMSPT')[p][0]=='M':
                pmatch_mloc[w] = p
                pmatches_m[w] = 1
            else:
                pmatch_fgkloc[w] = p
                pmatches_fgk[w] = 1
print 'Praesepe:',match_count

pmatch_loc = pmatch_loc[pmatch_loc>=0]
pmatch_mloc = pmatch_mloc[pmatch_mloc>=0]
pmatch_fgkloc = pmatch_fgkloc[pmatch_fgkloc>=0]


hlen=len(hdat)
hra = hdat.field('RA')
hdec = hdat.field('DEC')
heqw,hueqw = hdat.field('AVG_EQW'),hdat.field('AVG_EQW_ERR')
hll,hull = hdat.field('AVG_LHA'),hdat.field('AVG_LHA_ERR')
hperiods = hdat.field('PERIOD')
h_ppmax = calc_pmax.p_pmax(hdat.field('KH_MASS'),hperiods)
hbin = (hdat.field('BINARY')>0)
hpmem = hdat.field('ROESER_PMEM')

wh = at.read('../xray/wright_hyades_catalog.csv')
wra = wh['ra']
wdec = wh['dec']
wlen = len(wh)
hmatch_loc = np.ones(wlen,int)*-99
hmatch_mloc = np.ones(wlen,int)*-99
hmatch_fgkloc = np.ones(wlen,int)*-99
hmatches_m = np.zeros(wlen,int)
hmatches_fgk = np.zeros(wlen,int)

match_tol = 5./3500.
match_count = 0
for h in range(hlen):
    for w in range(wlen):
        if ((abs(wra[w]-hra[h])<match_tol) and (abs(wdec[w]-hdec[h])<match_tol)
            and ((hpmem[h]>=pmem_threshold) or (hpmem[h]<0)) 
            and (hbin[h]==False)
            ):
            #print 'found!',hdat.field('TWOMASSNAME')[h],
            #print hdat.field('KH_MASS')[h],hdat.field('MDM_SPECMATCH')[h]
            match_count += 1
            hmatch_loc[w] = h
            if hdat.field('MDM_SPEC_ADAMSPT')[h][0]=='M':
                hmatch_mloc[w] = h
                hmatches_m[w] = 1
            else:
                hmatch_fgkloc[w] = h
                hmatches_fgk[w] = 1
print 'Hyades:',match_count


#hmatch_loc = hmatch_loc[hmatch_loc>=0]
hmatch_mloc = hmatch_mloc[hmatch_mloc>=0]
hmatch_fgkloc = hmatch_fgkloc[hmatch_fgkloc>=0]

goodp = np.where(wprae['mass']<=0.8)[0]
goodh = np.where(wh['mass']<=0.8)[0]

# ChaMP/ChESS stars

allchamp = pyfits.open('../xray/ChaMPstars-final.fits')
cdat = allchamp[1].data

tab4 = at.read('../xray/tab4.dat')
spt = np.empty(len(tab4),'S1')
spt[:] = tab4['Type'][:]
eqws = -1*tab4['EqW']
num4 = len(tab4)
subtype = np.ones(num4,int)*-99
LxLbol = np.ones(num4)*-99.
LxLbolerr = np.ones(num4)*-99.
imJ = np.ones(num4)*-99.
imJerr = np.ones(num4)*-99.
LhaLbol_pub = 10**(tab4['LhaLbol'])
ulim_flag = np.zeros(num4,int)

tab2 = at.read('../xray/tab2.dat')

tab3 = at.read('../xray/tab3.dat')

for i in range(num4):
    j = np.where(tab2['CXOMP']==tab4['CXOMP'][i])[0]
    k = np.where(tab3['CXOMP']==tab4['CXOMP'][i])[0]
    subtype[i] = tab4['Type'][i][1]
    #print j, k
    if len(j)==1:
        LxLbol[i] = 10**(tab2['logXs/B'][j])
        LxLbolerr[i] = abs(tab2['e_logXs/B'][j]*np.log(10.)*LxLbol[i])
    else:
        print 'multi/no match!', j, tab4['CXOMP'][i]
    if len(k)==1:
        imJ[i] = tab3['imag'][k]-tab3['Jmag'][k]
        imJerr[i] = np.sqrt(tab3['e_imag'][k]**2 + tab3['e_Jmag'][k]**2)
    else:
        print 'multi/no match!', k, tab4['CXOMP'][i]
    m = np.where(cdat.field('XID')==('CXOMP'+tab4['CXOMP'][i]))[0]
    if len(m)==1:
        ulim_flag[i] = cdat.field('SPEC_HAUPPERLIMIT')[m]
    else:
        print 'multi/no match!', m, tab4['CXOMP'][i]

chi,err_chi = calc_chi.chi(imJ,'i-J')
LhaLbol = -1*eqws*(10**chi)
#LhaLbol[eqws<0] = -99.0
#LhaLbol[chi<0] = -99.0


plt.figure(figsize=(8,12))
ax = plt.subplot(211)
ax.plot(LxLbol[(spt=='M') & (ulim_flag==0)],eqws[(spt=='M') & (ulim_flag==0)],'s',
    mec='Gray',mfc='None',mew=2)
ax.plot(LxLbol[(spt!='M') & (ulim_flag==0)],eqws[(spt!='M') & (ulim_flag==0)],'+',
    mec='Gray',mew=2,ms=10)
ax.plot(10**(wprae['LxLbol'][pmatches_m==1]),-1*peqw[pmatch_mloc],'s',
    mec='DarkBlue',mfc='None',mew=2)
ax.plot(10**(wh['LxLbol'][hmatches_m==1]),-1*heqw[hmatch_mloc],'s',
    mec='#CF3800',mfc='None',mew=2)
ax.plot(10**(wprae['LxLbol'][pmatches_fgk==1]),-1*peqw[pmatch_fgkloc],'+',
    mec='DarkBlue',mew=2,ms=10)
ax.plot(10**(wh['LxLbol'][hmatches_fgk==1]),-1*heqw[hmatch_fgkloc],'+',
    mec='#CF3800',mew=2,ms=10)
ax.set_xscale('log')
ax.set_ylim(5,-16)
ax.set_ylabel(r'H$\alpha$ EqW',fontsize='x-large')
ax.tick_params(which='both',width=2,labelsize='x-large')
textx = 1.5e-6
ax.text(textx,-10,'Field',color='Gray',fontsize='large')
ax.text(textx,-8.75,'Praesepe',color='DarkBlue',fontsize='large')
ax.text(textx,-7.5,'Hyades',color='#CF3800',fontsize='large')
ax.plot(1,100,'ks',mew=2,ms=10,mfc='None',label='M dwarfs')
ax.plot(1,100,'k+',mew=2,ms=10,label='FGK dwarfs')
ax.legend(loc=2,numpoints=1,frameon=False)
ax.set_xlim(1e-6,1e-2)

#upper limits
up_lim = np.zeros(16).reshape(2,-1)
up_lim[0] = -1.1
ax.errorbar(LxLbol[(spt=='M') & (ulim_flag==1)],eqws[(spt=='M') & (ulim_flag==1)],up_lim,fmt=None,
    ecolor='DarkGrey',color='DarkGrey',lolims=True)




pll[((peqw-pueqw)<=0) & (peqw>0)] = -99.9

ax = plt.subplot(212)
#ax.plot(LxLbol[(spt=='M') & (ulim_flag==0)],LhaLbol[(spt=='M') & (ulim_flag==0)],'s',
#    mec='Gray',mfc='None',mew=2)
#ax.plot(LxLbol[(spt!='M') & (ulim_flag==0)],LhaLbol[(spt!='M') & (ulim_flag==0)],'+',
#    mec='Gray',mew=2,ms=10)
ax.errorbar(LxLbol[(spt=='M') & (ulim_flag==0)],
    LhaLbol[(spt=='M') & (ulim_flag==0)],
    xerr=LxLbolerr[(spt=='M') & (ulim_flag==0)],fmt='s',
    ecolor='Gray',mec='Gray',mfc='None',mew=2)
ax.errorbar(LxLbol[(spt!='M') & (ulim_flag==0)],
    LhaLbol[(spt!='M') & (ulim_flag==0)],
    xerr=LxLbolerr[(spt!='M') & (ulim_flag==0)],fmt='+',
    ecolor='Gray',mec='Gray',mew=2,ms=10)


ax.plot(10**(wprae['LxLbol'][pmatches_m==1]),pll[pmatch_mloc],'s',
    mec='DarkBlue',mfc='None',mew=2)
ax.plot(10**(wh['LxLbol'][hmatches_m==1]),hll[hmatch_mloc],'s',
    mec='#CF3800',mfc='None',mew=2)
ax.plot(10**(wprae['LxLbol'][pmatches_fgk==1]),pll[pmatch_fgkloc],'+',
    mec='DarkBlue',mfc='None',mew=2,ms=10)
ax.plot(10**(wh['LxLbol'][hmatches_fgk==1]),hll[hmatch_fgkloc],'+',
    mec='#CF3800',mfc='None',mew=2,ms=10)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylim(1e-6,1e-3)
ax.set_xlim(1e-6,1e-2)
ax.set_xlabel(r'$L_X/L_{bol}$',fontsize='xx-large')
ax.set_ylabel(r'$L_{H\alpha}/L_{bol}$',fontsize='xx-large')
ax.tick_params(which='both',width=2,labelsize='x-large')

#upper limits
up_lim = np.zeros(16).reshape(2,-1)
up_lim[0] = -1.1
ax.errorbar(LxLbol[(spt=='M') & (ulim_flag==1)],LhaLbol[(spt=='M') & (ulim_flag==1)],LhaLbol[(spt=='M') & (ulim_flag==1)]*0.3,fmt=None,
    ecolor='DarkGrey',color='DarkGrey',lolims=True)





plot_color='DarkBlue'
ulim_lha, uperr_ulim_lha, dnerr_ulim_lha = simple(peqw,pueqw,pdat.field('CHI'))
# No upper limits on M dwarfs for this set
good_stars = pmatch_fgkloc
pup = (np.ones(2*len(good_stars))*1e-6).reshape(2,-1)
pup[1] = uperr_ulim_lha[good_stars]
pup[0] = dnerr_ulim_lha[good_stars]
ax.errorbar(10**wprae['LxLbol'][pmatches_fgk==1],
    ulim_lha[good_stars],
    pup,ecolor=plot_color,
    color=plot_color,fmt=None,lolims=True)

plot_color='#CF3800'
ulim_lha, uperr_ulim_lha, dnerr_ulim_lha = simple(heqw,hueqw,hdat.field('CHI'))
# No upper limits on M dwarfs for this set
good_stars = hmatch_fgkloc
pup = (np.ones(2*len(good_stars))*1e-6).reshape(2,-1)
pup[1] = uperr_ulim_lha[good_stars]
pup[0] = dnerr_ulim_lha[good_stars]
ax.errorbar(10**wh['LxLbol'][hmatches_fgk==1],
    ulim_lha[good_stars],
    pup,ecolor=plot_color,
    color=plot_color,fmt=None,lolims=True)

plt.savefig('paperxray.png')
plt.savefig('paperxray.ps')
