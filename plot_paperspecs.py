# Script to plot sample spectra for our data
# Stephanie Douglas 11 June 2013

import pyfits, get_data, read_spec
from praesepe_comp import getspt
from makemodel import falt2
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np

pdat,pobs,pobsnr,pobsr = get_data.get_data('P')
hdat,hobs,hobsnr,hobsr = get_data.get_data('H')
pspts = getspt(pdat.field('ADAMSPT'))
hspts = getspt(hdat.field('MDM_SPEC_ADAMSPT'))

resols = {'mdm':6.,'hydra':3.3,'mage':0.8}

def plot_spec(ax,filename,source,offset,style='k-',smooth=True):
    w,f,n = read_spec.get_spec(filename,source)
    if smooth:
        f = falt2(w,f,resols[source]*0.75)
    f_norm = f/np.average(f[(w>=6548) & (w<=6558)])
    ax.plot(w,f_norm+offset,style,lw=0.75)


# Hyades - MDM - M2
hyad_act = 'Praesepe_Hyades/MODspec/Nov_2012_specs/2012_11_11_finals_steph/04360416+185318.fits'
hyad_inact = 'Praesepe_Hyades/MODspec/Nov_2012_specs/2012_11_11_finals_steph/04363893+183656.fits'#

# MDM - M3 for Praesepe
mdm_act = 'Praesepe_Hyades/MODspec/Feb2011_n4/HSHJ303.fits'
mdm_inact = 'Praesepe_Hyades/MODspec/Feb2011_n4/HSHJ272.fits'

# Hydra
#hyd_act = 'Hydra_spectra/PraeB1_ascii/JS430-273.asc'
#hyd_inact = 'Hydra_spectra/PraeF1_ascii/JS365-509.asc'
hyd_act = 'Hydra_spectra/PraeB1_new/JS430-273.fits'
hyd_inact = 'Hydra_spectra/PraeF1_new/JS365-509.fits'

# Mage 
mage_act = 'Praesepe_Hyades/MagE/JS123_F.fits'



textx = 5000
plt.figure(figsize=(12,6))
ax4 = plt.subplot2grid((1,5),(0,0),colspan=4)
plot_spec(ax4,hyad_inact,'mdm',6,'r-')
plot_spec(ax4,hyad_act,'mdm',6)
ax4.text(textx,7.55,'ModSpec (Hyades)',fontstyle='oblique')
ax4.text(textx,7.25,'2MASS J0436+1853')
ax4.text(textx,6.95,'2MASS J0436+1836',color='r')
plot_spec(ax4,mdm_inact,'mdm',3.9,'r-')
plot_spec(ax4,mdm_act,'mdm',4)
ax4.text(textx,5.25,'ModSpec (Praesepe)',fontstyle='oblique')
ax4.text(textx,4.95,'JS506')
ax4.text(textx,4.65,'JS447',color='r')
plot_spec(ax4,hyd_inact,'hydra',2,'r-')
plot_spec(ax4,hyd_act,'hydra',2)
ax4.text(textx,3.25,'Hydra (Praesepe)',fontstyle='oblique')
ax4.text(textx,2.95,'JS430')
ax4.text(textx,2.65,'JS365',color='r')
plot_spec(ax4,mage_act,'mage',0)
ax4.text(textx,1.25,'MagE (Praesepe)',fontstyle='oblique')
ax4.text(textx,0.95,'JS298')#'JS123')
ax4.set_xlim(4500,8250)
ax4.set_ylim(0,8.25)
ax4.set_xlabel('Wavelength ($\AA$)',fontsize='x-large')
ax4.set_ylabel('flux / flux(6555 $\AA$) + offset',fontsize='x-large')
ax4.tick_params(labelsize='large')

ax5 = plt.subplot2grid((1,5),(0,4))
ax5.add_patch(Rectangle((6548,0),10,11,fc='#DCDCDC',ec='none',fill=True))
ax5.add_patch(Rectangle((6570,0),10,11,fc='#DCDCDC',ec='none',fill=True))
plot_spec(ax5,hyad_inact,'mdm',6,'r-',False)
plot_spec(ax5,hyad_act,'mdm',6,smooth=False)
plot_spec(ax5,mdm_inact,'mdm',3.9,'r-',False)
plot_spec(ax5,mdm_act,'mdm',4,smooth=False)
plot_spec(ax5,hyd_inact,'hydra',2,'r-',False)
plot_spec(ax5,hyd_act,'hydra',2,smooth=False)
plot_spec(ax5,mage_act,'mage',0,smooth=False)
ax5.set_xlim(6475,6625)
ax5.set_ylim(0,8.25)
ax5.set_xticks((6500,6600))
ax5.tick_params(labelsize='large',labelleft='off')
plt.savefig('paper_specex.png')
plt.savefig('paper_specex.ps',orientation='landscape')
#savefig('paper_specha.png')
#savefig('paper_specha.ps')


"""
M0 = where((hspts>=20) & (hspts<21.))[0]
M1 = where((pspts>=21) & (pspts<22.) & (pdat.field('MDM_SPECMATCH')>0))[0]
M2 = where((pspts>=22) & (pspts<23.) & (pdat.field('HYDRA_SPECMATCH')>0))[0]
M3 = where((pspts>=23) & (pspts<24.) & (pdat.field('MAGE_SPECMATCH')>0))[0]

offset=0

for i in M0:
    w,f,n = read_spec.get_spec(hdat.field('MDM_SPEC_FILE')[i],'mdm',basepath='/home/stephanie/ptf/spectra/Praesepe_Hyades/MODspec/')
    f_norm = f/average(f[(w>=6548) & (w<=6558)])
    plot(w,f_norm+offset,'k-')
    offset+=1
    print hdat.field('MDM_SPEC_FILE')[i],i

for i in M1:
    w,f,n = read_spec.get_spec(pdat.field('MDM_SPEC_NAME')[i],'mdm')
    f_norm = f/average(f[(w>=6548) & (w<=6558)])
    plot(w,f_norm+offset,'k-')
    offset+=1
    print pdat.field('MDM_SPEC_NAME')[i], i

for i in M2:
    w,f,n = read_spec.get_spec(pdat.field('HYDRA_FILE')[i],'hydra')
    f_norm = f/average(f[(w>=6548) & (w<=6558)])
    plot(w,f_norm+offset,'k-')       
    offset+=1
    print pdat.field('HYDRA_FILE')[i], i

for i in M3:
    w,f,n = read_spec.get_spec(pdat.field('MAGE_SPEC_FILE')[i],'mage')    f_norm = f/average(f[(w>=6548) & (w<=6558)])
    plot(w,f_norm+offset,'k-')       
    offset+=2
    print pdat.field('MAGE_SPEC_FILE')[i], i
"""

