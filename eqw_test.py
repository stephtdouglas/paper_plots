import numpy as np
import matplotlib.pyplot as plt
import get_data
import asciitable as at

#for Praesepe
pdat,pobs,pobsnr,pobsr = get_data.get_data('P')
kh = at.read('../../catalogs/KafkaHoneycutt_halpha.tsv')
all_kh_eqws = kh['EW(Ha)']

# going to re-write because I think the matching here is counter-intertuitive. Want to require that a star has a measurement in our catalog before seeking out value in K&H
our_kh = np.where(pdat.field('KAFKAMATCH')>0)[0]
our_halpha = pdat.field('KAFKA_EQW')[our_kh]
our_halpha_err = pdat.field('KAFKA_EQW_ERR')[our_kh]
our_color = pdat.field('RPRIME_K')[our_kh]
our_ra = pdat.field('RA')[our_kh]
our_dec = pdat.field('DEC')[our_kh]

separation = np.ones(len(our_kh))*-99.
their_halpha = np.ones(len(our_kh))*-99.

match_tol = 5.
for i in range(len(our_kh)):
	for j in range(len(all_kh_eqws)):
		separation[i] = np.sqrt(((our_ra[i]-kh['_RA'][j])*np.cos(our_dec[i]*np.pi/180.)*3600.)**2+((our_dec[i]-kh['_DE'][j])*3600.)**2)
		if (separation[i]<match_tol):
			their_halpha[i] = all_kh_eqws[j]

# plot the straightforward comparison of EqW measurements

plt.figure()
## plt.plot(our_halpha*-1,kh['EW(Ha)'],'ko')
y = np.arange(-10,2,0.001)
x = np.arange(-10,2)
plt.plot(y, y, color='k')
##plt.Line2D(x, y, linewidth='2', ls='-', color = 'b')
lines = {'linestyle': 'None'}
yerr = our_halpha_err
xerr = 0.25
##plt.errorbar(our_halpha*-1, our_halpha_err, marker='o')
##plt.errorbar(our_halpha*-1, kh['EW(Ha)'], yerr, xerr, capsize=0, ls='none', color='red', marker = 'o') #elinewidth=2) 
plt.errorbar(their_halpha, our_halpha*-1,yerr,  xerr,  capsize=0, ls='none', color='DarkBlue', marker = 'o', markeredgecolor='DarkBlue') #elinewidth=2) 
plt.ylabel('Measured EqW ($\AA$)', fontsize='x-large')
plt.xlabel('Kafka & Honeycutt 2006 EqW ($\AA$)', fontsize='x-large')
plt.ylim(1.,-10.)
plt.xlim(1.,-10.)
ax = gca()
ax.tick_params(labelsize='large')
plt.show()
#plt.text(-2,-8,'161 stars')
#plt.title('Praesepe')
plt.legend(loc=2,numpoints=1,scatterpoints=1,frameon=False)
plt.savefig('eqw_praesepe_v2.eps', bbox_inches='tight')


#for Hyades
pdat_h,pobs,pobsnr,pobsr = get_data.get_data('H')
hya = at.read('MasterHyadesList.tsv')
stauffer91_eqws = hya['S91']
stauffer94_eqws = hya['S94']
stauffer97_eqws = hya['S97']
tendrup00_eqws = hya['T00']
jones96_eqws = hya['J96']

# let's try the same approach
our_h_ra = pdat_h.field('RA')
our_h_dec = pdat_h.field('Dec')
# our MDM data have emission > 0 A, which is bad. Need to *-1.
our_h_halpha = pdat_h.field('MDM_EQW')*-1.
our_h_halpha_err = pdat_h.field('MDM_EQW_ERR')
our_h_color = pdat_h.field('RPRIME_K')

separation_h = np.ones(len(our_h_ra))*-99.
stauffer91_matches = np.ones(len(our_h_ra))*-9999.

for i in range(len(our_h_ra)):
	for j in range(len(stauffer91_eqws)):
		separation_h[i] = np.sqrt(((our_h_ra[i]-hya['RA'][j])*np.cos(our_h_dec[i]*np.pi/180.)*3600.)**2+((our_h_dec[i]-hya['Dec'][j])*3600.)**2)
		if ((separation_h[i]<match_tol) and our_h_halpha[i]!=99. and stauffer91_eqws[j]!=-9999.):
			# Stauffer et al. 1991 also have emission > 0 A, so need to *-1.
			stauffer91_matches[i] = stauffer91_eqws[j]*-1.

# can apply same approach to other data
stauffer94_matches = np.ones(len(our_h_ra))*-9999.
stauffer97_matches = np.ones(len(our_h_ra))*-9999.
tendrup00_matches = np.ones(len(our_h_ra))*-9999.

for i in range(len(our_h_ra)):
	for j in range(len(stauffer97_eqws)):
		separation_h[i] = np.sqrt(((our_h_ra[i]-hya['RA'][j])*np.cos(our_h_dec[i]*np.pi/180.)*3600.)**2+((our_h_dec[i]-hya['Dec'][j])*3600.)**2)
		if ((separation_h[i]<match_tol) and our_h_halpha[i]!=99. and stauffer97_eqws[j]!=-9999.):
			# Stauffer et al. 1997 also have emission > 0 A, so need to *-1.
			stauffer97_matches[i] = stauffer97_eqws[j]*-1.

# there's only one match to Stauffer 1994... 
for i in range(len(our_h_ra)):
	for j in range(len(stauffer94_eqws)):
		separation_h[i] = np.sqrt(((our_h_ra[i]-hya['RA'][j])*np.cos(our_h_dec[i]*np.pi/180.)*3600.)**2+((our_h_dec[i]-hya['Dec'][j])*3600.)**2)
		if ((separation_h[i]<match_tol) and our_h_halpha[i]!=99. and stauffer94_eqws[j]!=-9999.):
			# Stauffer et al. 1994 also have emission > 0 A, so need to *-1.
			stauffer94_matches[i] = stauffer94_eqws[j]

match = np.where(stauffer94_matches!=-9999.)[0]


# plot the straightforward comparison of EqW measurements
good = np.where(our_h_halpha!=-99.)[0]
nspec = len(np.where(stauffer91_matches!=9999.)[0])+len(np.where(stauffer97_matches!=9999.)[0])+1 #that one is the Stauffer 1994 star
plt.figure()
x = np.arange(-10,2,0.001)
plt.plot(x, x, color='k')
lines = {'linestyle': 'None'}
yerr = our_h_halpha_err
xerr = 0.026  # quoted 1 sigma uncertainty in Stauffer 1991
plt.errorbar(stauffer91_matches[good], our_h_halpha[good],yerr[good], xerr,  capsize=0, ls='none', color='OrangeRed', marker = 'D', mec="OrangeRed", label='Stauffer+ 1991') #elinewidth=2) 
plt.errorbar(stauffer94_matches[match]*-1., our_h_halpha[match],yerr[match], xerr,  capsize=0, ls='none', color='OrangeRed', marker = '^', mec="OrangeRed", label='Stauffer+ 1994') #elinewidth=2) 
#print our_h_halpha[match], xerr[match]
plt.errorbar(stauffer97_matches[good], our_h_halpha[good], yerr[good],xerr,  capsize=0, ls='none', color='OrangeRed', marker = 's', mec="OrangeRed", label='Stauffer+ 1997') #elinewidth=2) 
plt.ylabel(r'Measured H$\alpha$ EqW ($\AA$)',fontsize='x-large')
plt.xlabel(r'Literature H$\alpha$ EqW ($\AA$)',fontsize='x-large')
ax=plt.gca()
ax.tick_params(labelsize='large')
plt.ylim(1.5,-8.75)
plt.xlim(1.5,-8.75)
#print nspec
ax.legend(numpoints=1,handletextpad=0.2,handlelength=1,borderaxespad=0.2,loc='upper left')
#plt.title('Hyades')
plt.savefig('eqw_Hyades.eps', bbox_inches='tight')
#plt.show()

