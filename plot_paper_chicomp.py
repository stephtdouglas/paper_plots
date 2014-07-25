# Script to calculate Chi=f_halpha/f_bol as a function of r-K color
# For a subset of the PMSU spectroscopic survey
#
# Created by Stephanie Douglas, 21 February 2013
# Updated/rewritten 23 April 2013
# and again 7 June and 11 June 2013
# and 17 july and 15 august
################################################################################

import numpy as np
import asciitable as at
from ha_cont import ha_cont
import os, bol_corr, read_spec, get_data, pickle
import scipy.odr as odr
#from scipy.optimize import minimize
import pyfits
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from praesepe_comp import getspt
from emissionline import emissionline

infil = open('/home/stephanie/ptf/chi/SDSSmatches.pkl','rb')
cmatch = pickle.load(infil)
infil.close()

sdssm = pyfits.open('/home/stephanie/ptf/chi/DR7_tot_4_0.fits')
plate = sdssm[1].data.field('PLATE')
mjd = sdssm[1].data.field('MJD')
fiber = sdssm[1].data.field('FIBER')

covey = pyfits.open('/home/stephanie/ptf/chi/superclean.fits')

cmatches = cmatch[cmatch>=0] # index in Andrew's file
cmatch_loc = np.where(cmatch>=0)[0] # Index in Kevin's file
cnum = len(cmatches)

m_cont = np.zeros(cnum)
m_fbol = np.zeros(cnum)
m_bcK = np.zeros(cnum)
m_eqw = np.ones(cnum)*-99.
m_cont2 = np.zeros(cnum)
m_cont3 = np.zeros(cnum)

m_r = covey[1].data.field('RMAG')[cmatch_loc]
m_rmK = m_r - covey[1].data.field('KMAG')[cmatch_loc]
m_rmKerr = np.sqrt(covey[1].data.field('RERR')**2+
     covey[1].data.field('KERR')**2)[cmatch_loc]
m_Kerr = covey[1].data.field('KERR')[cmatch_loc]

for i in range(cnum):
    j = cmatches[i]
    # Read in the SDSS spectrum
    if fiber[j] <10:
        fib = '00'+str(fiber[j])
    elif fiber[j]<100:
        fib = '0'+str(fiber[j])
    else:
        fib = str(fiber[j])
    if plate[j] < 10:
        plat = '000'+str(plate[j])
    elif plate[j]<100:
        plat = '00'+str(plate[j])
    elif plate[j]<1000:
        plat = '0'+str(plate[j])
    else:
        plat = str(plate[j])
    sdssfile = 'SDSSM/spSpec-%d-%s-%s.fit' % (mjd[j],plat,fib)
    w,f,n = read_spec.get_spec(sdssfile,'sdss')

    # Measure f_6560
    m_cont[i] = ha_cont(w,f,start_cont=6550,end_cont=6560
        ,start_cont2=6572,end_cont2=6582) # erg/s/cm2

    m_eqw[i],junk1,junk2,junk3 = emissionline(w,f,n,
#        6548,6558,6570,6580,6558,6570)
        6550,6560,6572,6582,6560,6572)

    # Calculate fbol
    m_fbol[i] = bol_corr.kh_fbol(m_r[i],
        m_rmK[i],10,'r','r','K')
    m_bcK[i] = bol_corr.kraus(m_rmK[i])

    # Testing other continuum regions
    m_cont2[i] = ha_cont(w,f,start_cont=6555,end_cont=6560) # erg/s/cm2
    m_cont3[i] = ha_cont(w,f,start_cont=6550,end_cont=6560
        ,start_cont2=6570,end_cont2=6580) # erg/s/cm2

m_chi2 = np.log10(m_cont2/m_fbol)
m_chi3 = np.log10(m_cont3/m_fbol)



m_eqw = m_eqw*-1.
m_chi1 = m_cont/m_fbol
m_err_part = m_rmKerr**2 * m_bcK**2 / m_rmK**2 + m_Kerr**2
m_err_chi1 = m_chi1*(0.4*np.log(10.))*np.sqrt(m_err_part) #chi*0.04
m_err_chi = np.log10(np.e)*m_err_chi1/m_chi1
m_chi = np.log10(m_chi1)

andrew_chi = np.log10(sdssm[1].data.field('LHALBOL')/
    sdssm[1].data.field('EWHA'))[cmatches]
m_spt = sdssm[1].data.field('SPT')[cmatches]+48
#m_chi[isnan(m_chi)] = 99.0
#m_err_chi[isnan(m_chi)] = 50.0

bad = np.where(np.isnan(m_chi))[0]
m_chi = np.delete(m_chi,bad)
m_err_chi = np.delete(m_err_chi,bad)



def linear(p,x):
    a,b=p
    return a*x + b

def quad(p,x):
    a,b,c = p
    return a*x**2 + b*x + c

def cube(p,x):
    a,b,c,d = p
    return a*x**3 + b*x**2 + c*x + d

def quart(p,x):
    a,b,c,d,g = p
    return a*x**4 + b*x**3 + c*x**2 + d*x + g

def fitit(order,xdata,ydata,xerrs,yerrs,verb=True):
    max_flag = 0
    if order==1:
        mymod = odr.Model(linear)
        maxit = 200
    elif order==2:
        mymod = odr.Model(quad)
        maxit = 400
    elif order==3:
        mymod = odr.Model(cube)
        maxit = 500
    elif order==4:
        mymod = odr.Model(quart)
        maxit = 800
    p0 = np.ones(order+1)
    mydata = odr.RealData(xdata,ydata,xerrs,yerrs)
    myodr = odr.ODR(mydata,mymod,beta0=p0,maxit=200)
    myoutput = myodr.run()
    if verb:
        myoutput.pprint()
    if myoutput.stopreason[0]=='Iteration limit reached':
        print myoutput.stopreason
        max_flag = 1
    return myoutput.beta,myoutput.sd_beta,max_flag


phx_cols = {'u':'SDSS_U','g':'SDSS_G','r':'SDSS_R','i':'SDSS_I','z':'SDSS_Z',
    'J':'TWO_J','K':'TWO_K'}

covey_cols = {'u':'UMAG','g':'GMAG','r':'RMAG','i':'IMAG','z':'ZMAG',
    'J':'JMAG','K':'KMAG'}

covey_errs = {'u':'UERR','g':'GERR','r':'RERR','i':'IERR','z':'ZERR',
    'J':'JERR','K':'KERR'}

prae_cols = {'u':'SDSSPSFMAG_U','g':'SDSSPSFMAG_G','r':'SDSSPSFMAG_G',
    'i':'SDSSPSFMAG_I','z':'SDSSPSFMAG_Z','J':'TWOMASS_J','K':'TWOMASS_K'}

prae_errs = {'u':'SDSSPSFMAGERR_U','g':'SDSSPSFMAGERR_G','r':'SDSSPSFMAGERR_G',
    'i':'SDSSPSFMAGERR_I','z':'SDSSPSFMAGERR_Z','J':'TWOMASS_JERR',
    'K':'TWOMASS_KERR'}
ckeys = ['u','g','r','i','z','J','K']
carr = np.array(ckeys)
cdat = covey[1].data


m_spt = np.delete(m_spt,bad)
m_chi = 10**m_chi
m_chi2 = np.delete(m_chi2,bad)
m_chi2 = 10**m_chi2


avg_chi = np.zeros(10)
std_chi = np.zeros(10)
avg_chi2 = np.zeros(10)
std_chi2 = np.zeros(10)
for i in np.arange(48,58,1):
    good = np.where(m_spt==i)[0]
    avg_chi[i-48] = np.average(m_chi[good])
    std_chi[i-48] = np.std(m_chi[good])
    avg_chi2[i-48] = np.average(m_chi2[good])
    std_chi2[i-48] = np.std(m_chi2[good])

#plt.errorbar(arange(48,58),avg_chi,std_chi,fmt='o',
#    label='Calculations for SDSS stars',color='k',ms=8,
#    elinewidth=2)
plt.errorbar(np.arange(48,58),avg_chi2,std_chi2,fmt='o',
    label='Calculations for SDSS stars',color='k',ms=8,
    elinewidth=2,capsize=0)

w_chi2 = np.array([1.16e-4,1.16e-4,0.966e-4,0.738e-4,0.637e-4,0.274e-4,0.176e-4,
    0.052e-4,0.06e-4,0.038e-4])
w_chi2_err = np.array([0.277e-4,0.451e-4,0.306e-4,0.216e-4,0.286e-4,0.128e-4,
    0.052e-4,0.015e-4,0.029e-4,0.011e-4])
plt.errorbar(np.arange(48,58)+0.05,w_chi2,w_chi2_err,fmt='o',
    label='WHW04 values given in West&Hawley08',color='g',mew=2,
    ms=8,elinewidth=2,mec='g',mfc='None',capsize=0)

#compare with spectral-type relations from Schmidt et al. 2014 (ArXiv 1406.1228v1)
#s_spt = np.arange(7,10,0.5)
#s14_spt = s_spt+(55-7)
#s14_chi = 1.09*np.exp(-1.0*s_spt/0.602) - (1.74e-7)*s_spt + 3.88e-6
#s14_chi1 = s14_chi + 1.43e-6
#s14_chi2 = s14_chi - 1.43e-6

ax = plt.gca()
#ax.plot(s14_spt,s14_chi,'r-.',lw=3,label='Schmidt+14')
#ax.plot(s14_spt,s14_chi1,'r-.',lw=1.5)
#ax.plot(s14_spt,s14_chi2,'r-.',lw=1.5)

plt.yscale('log')
plt.legend(loc=3,numpoints=1)
plt.xlim(47.5,57.5)
plt.ylim(1e-6,2e-4)
ax.set_xticks(np.arange(48,58))
ax.set_xticklabels(['M0','M1','M2','M3','M4','M5','M6','M7','M8','M9'])
ax.tick_params(labelsize='large')
ax.set_xlabel('Spectral Type',fontsize='large')
ax.set_ylabel(r'$\chi$',fontsize='x-large')


plt.savefig('paper_chicomp.png')
plt.savefig('paper_chicomp.eps',bbox_inches='tight')
