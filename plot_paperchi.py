# Script to calculate Chi=f_halpha/f_bol as a function of r-K color
# For a subset of the PMSU spectroscopic survey
#
# Created by Stephanie Douglas, 21 February 2013
# Updated/rewritten 23 April 2013
# and again 7 June and 11 June 2013
################################################################################

import numpy as np
import asciitable as at
from ha_cont import ha_cont
import os, bol_corr, read_spec, get_data, pickle, calc_chi
import scipy.odr as odr
#from scipy.optimize import minimize
import pyfits
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from praesepe_comp import getspt
from emissionline import emissionline



# Get matched photometry data
phot = at.read('../../spectra/PMSU/pmsu_phot.tsv')
# (r-K)
rmK = phot['r']-phot['K']
add_err = 0.1
rmKerr = np.sqrt((phot['r_err']/100.0)**2 + phot['K_err']**2 + add_err**2)
Kerr = phot['K_err']**2

pnum = len(phot)

# Empty arrays for data
hacont = np.zeros(pnum)
fbol = np.zeros(pnum)
bcK = np.zeros(pnum)
eqw = np.ones(pnum)*-99.0

#figure()
# Loop over everything to calculate fluxes
for i in range(pnum):
    # Read in the PMSU spectrum
    pmsufile = ('/home/stephanie/ptf/spectra/PMSU/book%d-1.txt' %
                        phot['pmsu'][i])
    w,f,n = read_spec.get_spec(pmsufile,'pmsu')

    # Measure f_6560
    hacont[i] = ha_cont(w,f,start_cont=6548,end_cont=6558
        ,start_cont2=6570,end_cont2=6580) # erg/s/cm2

    eqw[i],junk1,junk2,junk3 = emissionline(w,f,n,
        6548,6558,6570,6580,6558,6570)

    # Calculate fbol
    fbol[i] = bol_corr.kh_fbol(phot['r'][i],
        rmK[i],10,'r','r','K')
    bcK[i] = bol_corr.kraus(rmK[i])

badi = np.where(rmKerr>1)[0]
rmK = np.delete(rmK,badi)
rmKerr = np.delete(rmKerr,badi)
Kerr = np.delete(Kerr,badi)
hacont = np.delete(hacont,badi)
fbol = np.delete(fbol,badi)
bcK = np.delete(bcK,badi)
eqw = np.delete(eqw,badi)
eqw = eqw*-1.

# Calculate chi
chi1 = hacont/fbol
err_part = rmKerr**2 * bcK**2 / rmK**2 + Kerr**2
err_chi1 = chi1*(0.4*np.log(10.))*np.sqrt(err_part) #chi*0.04
err_chi = np.log10(np.e)*err_chi1/chi1
chi = np.log10(chi1)

# Read in my chi calculations from the PHOENIX grid
phx = pyfits.open('../../chi/Phoenix_Steph.fits')
p_chi = phx[1].data.field('CHI')
p_rpmK = phx[1].data.field('RPK')
pd = phx[1].data
p_r = pd.field('SDSS_RP')
p_rmi = (pd.field('SDSS_R')-pd.field('SDSS_I'))
p_fbol = pd.field('FBOL')
p_hac = pd.field('HA_CONT')
p_eqw = pd.field('HA_EQW')
p_bcK = pd.field('BC_K')
p_err_chi = (p_hac/p_fbol)*(0.4*np.log(10.))*np.sqrt(0.02)

# Read in matches between Kevin's superclean catalog and Andrew's SDSS M dwarfs

infil = open('../../chi/SDSSmatches.pkl','rb')
cmatch = pickle.load(infil)
infil.close()

sdssm = pyfits.open('../../chi/DR7_tot_4_0.fits')
plate = sdssm[1].data.field('PLATE')
mjd = sdssm[1].data.field('MJD')
fiber = sdssm[1].data.field('FIBER')

covey = pyfits.open('../../chi/superclean.fits')

cmatches = cmatch[cmatch>=0] # index in Andrew's file
cmatch_loc = np.where(cmatch>=0)[0] # Index in Kevin's file
cnum = len(cmatches)

m_cont = np.zeros(cnum)
m_fbol = np.zeros(cnum)
m_bcK = np.zeros(cnum)
m_eqw = np.ones(cnum)*-99.

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

def linear(p,x):
    a,b=p
    return a*x + b

def quad(p,x):
    a,b,c = p
    return a*x**2 + b*x + c

def cube(p,x):
    a,b,c,d = p
    return a*x**3 + b*x**2 + c*x + d

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


ckeys = ['u','g','r','i','z','J','K']
carr = np.array(ckeys)
cdat = covey[1].data

def color(c1,c2):
    j = np.where(carr==c1)[0]
    k = np.where(carr==c2)[0]
    print c1,c2
    if k<j:
       junk = c1
       c1 = c2
       c2 = junk
       print c1,c2
    p_color = (pd.field(phx_cols[c1]) - pd.field(phx_cols[c2]))
    m_color = (cdat.field(covey_cols[c1]) - cdat.field(covey_cols[c2]))[cmatch_loc]
    m_color = np.delete(m_color,bad)
    m_error = np.sqrt(cdat.field(covey_errs[c1])**2 + cdat.field(covey_errs[c2])**2)[cmatch_loc]
    m_error = np.delete(m_error,bad)
    str_title = '%s-%s' % (c1,c2)
    return p_color,  m_color, m_error, str_title

m_eqw = m_eqw*-1.
m_chi1 = m_cont/m_fbol
m_err_part = m_rmKerr**2 * m_bcK**2 / m_rmK**2 + m_Kerr**2
m_err_chi1 = m_chi1*(0.4*np.log(10.))*np.sqrt(m_err_part) #chi*0.04
m_err_chi = np.log10(np.e)*m_err_chi1/m_chi1
m_chi = np.log10(m_chi1)
bad = np.where(np.isnan(m_chi))[0]
m_chi = np.delete(m_chi,bad)
m_err_chi = np.delete(m_err_chi,bad)


m_rmi = m_r - covey[1].data.field('IMAG')[cmatch_loc]
m_rmierr = np.sqrt(covey[1].data.field('RERR')[cmatch_loc]**2 + covey[1].data.field('IERR')[cmatch_loc]**2)
m_rp = m_r - 0.035*m_rmi - 0.007
m_rpmK = m_rp - covey[1].data.field('KMAG')[cmatch_loc]
m_rpmKerr = np.sqrt(covey[1].data.field('RERR')[cmatch_loc]**2 + (0.035*m_rmierr)**2 + covey[1].data.field('KERR')[cmatch_loc]**2)
m_rp = np.delete(m_rp,bad)
m_rpmK = np.delete(m_rpmK,bad)
m_rpmKerr = np.delete(m_rpmKerr,bad)

plt.figure()
ax = plt.subplot(111)
plt.plot(m_rpmK,m_chi,'s',label='SDSS',mfc='DarkGrey',mec='DarkGrey')
plt.plot(rmK,chi,'ko',label='PMSU')
plt.plot(p_rpmK[2:],p_chi[2:],'mx',label='PHOENIX',ms=8,mew=3)
x = np.arange(2,9,0.1)
plt.plot(x,calc_chi.chi(x)[0],'m-')
plt.legend(loc='best',numpoints=1)
plt.xlabel('(r\'-K)',fontsize='x-large')
plt.ylabel('log($\chi$)',fontsize='x-large')
plt.xlim(2,8.5)
plt.ylim(-6,-3.5)
ax.tick_params(labelsize='large')
plt.savefig('paperchi.png')
plt.savefig('paperchi.ps')



p_imJ, m_imJ, m_imJerr, st = color('i','J')
p_imJerr = p_imJ*0.05


w08_imJ1 = np.arange(1.61,3.22,0.001)
w08_imJ2 = np.arange(3.22,4.29,0.001)
w08_chi1 = (5.25e-4)*np.exp(-1.0*w08_imJ1/(1.53)) - 6.09e-5
w08_chi2 = (2.257e-3)*np.exp(-1.0*w08_imJ2/(2.01e3)) - 2.248e-3

# Matching first panel in Figure 2 in West & Hawley (2008)
plt.figure()
ax = plt.subplot(111)
ax.plot(m_imJ,m_chi,'s',label='SDSS',mfc='DarkGrey',mec='DarkGrey')
#p,up,fl = fitit(2,m_imJ,m_chi,m_imJerr,m_err_chi)
x = np.arange(1,4.5,0.5)
#ax.plot(x,quad(p,x),'b-',label='SDSS',lw=2)
#fit_str = '{0:.4f}(i-J)^2 + {1:.4f}(i-J) + {2:.4f}'.format(p[0],p[1],p[2])
#ax.text(1.1,-5.5,fit_str,color='b')
ax.plot(p_imJ[2:],p_chi[2:],'mx',label='PHOENIX',ms=8,mew=3)
#p,up,fl = fitit(2,p_imJ[2:],p_chi[2:],p_imJerr[2:],p_err_chi[2:])
#fit_str = '{0:.4f}(i-J)^2 + {1:.4f}(i-J) + {2:.4f}'.format(p[0],p[1],p[2])
ax.plot(x,calc_chi.chi(x,'i-J')[0],'m-',lw=2)

#calced_chi = calc_chi.chi(x,'i-J')[0]
#ax.plot(x,calced_chi + calced_chi*0.02,'m-',lw=1)
#ax.plot(x,calced_chi - calced_chi*0.02,'m-',lw=1)

#ax.text(1.1,-5.75,fit_str,color='m')
ax.set_ylabel('log($\chi$)',fontsize='x-large')
ax.set_xlabel("(i-J)",fontsize='x-large')

ax.plot(w08_imJ1,np.log10(w08_chi1),'g--',lw=2)
ax.plot(w08_imJ2,np.log10(w08_chi2),'g--',lw=2,label='West&Hawley08')

#compare the chi values given in Schmidt et al. 2014 (ArXiv 1406.1228v1)
#s14_imJ = np.arange(2.6,4.5,0.01)
#s14_chi = 1.39e-3 * np.exp(-1.0*s14_imJ/0.602) + 0.939e-6
#s14_chi1 = 1.39e-3 * np.exp(-1.0*s14_imJ/0.602) + 0.939e-6 + 0.499e-6
#s14_chi2 = 1.39e-3 * np.exp(-1.0*s14_imJ/0.602) + 0.939e-6 - 0.499e-6
#ax.plot(s14_imJ,np.log10(s14_chi),'r-.',lw=3,label='Schmidt+14')
#ax.plot(s14_imJ,np.log10(s14_chi1),'r-.',lw=1.5)
#ax.plot(s14_imJ,np.log10(s14_chi2),'r-.',lw=1.5)

plt.legend(loc='best',numpoints=1)
plt.savefig('paperchi_imJ.png')
plt.savefig('paperchi_imJ.ps')
