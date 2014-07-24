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

# now want to look at median difference in EqW as a function of color
#color_bins = np.arange(2.0,6.5,0.5)
#print color_bins
#binned_ha = np.zeros(8)
#binned_errors = np.zeros(8)
#diff_ha = (our_halpha*-1)-their_halpha
#diff_err = np.sqrt(our_halpha_err**2 + 0.25**2)

#for i in range(len(color_bins)-1):
#	in_this_bin = np.where((our_color>color_bins[i]) & 
#			(our_color<=color_bins[i+1]))[0]
#	binned_ha[i] = np.median(diff_ha[in_this_bin])
#	binned_errors[i] = np.median(diff_err[in_this_bin])

#print binned_ha, binned_errors

#fig = plt.figure()
#ax = fig.add_subplot(1,1,1)
##plt.plot(kh['EW(Ha)'],((our_halpha*-1)-kh['EW(Ha)']),'ro')
#plt.plot(our_color, diff_ha, 'bo')
##plt.errorbar(our_color,((our_halpha*-1)-kh['EW(Ha)']), our_halpha_err, capsize=0, ls='none', color='red', marker = 'o') #elinewidth=2) 
##ax.set_xscale('log')
#x = np.arange(-2,7)
#y = x/x - 1
#plt.plot(x, y, linewidth = 3, linestyle = '-')
#plt.axhline(y=0, color='green')
#plt.step(color_bins[:-1], binned_ha, where='post')
#plt.errorbar(color_bins[:-1]+0.25, binned_ha, binned_errors, capsize = 0, ls='none', color ='red', marker = '*', ms=15)
##plt.xlabel('K&H06 EqW')
#plt.xlabel('($r^{\prime}-K$)')
#plt.ylabel('Our EqW - K&H06 EqW ($\AA$)')
##plt.xlim(-8.5,0.5)
#plt.xlim(2,6.0)
#plt.ylim(-1.75,0.5)
#plt.title('Praesepe')
##plt.show()
#plt.savefig('eqw_praesepe_diff.eps', bbox_inches='tight')

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

# looks like there are no matches to the Tendrup 2000 data
#for i in range(len(our_h_ra)):
#	for j in range(len(tendrup00_eqws)):
#		separation_h[i] = np.sqrt(((our_h_ra[i]-hya['RA'][j])*np.cos(our_h_dec[i]*np.pi/180.)*3600.)**2+((our_h_dec[i]-hya['Dec'][j])*3600.)**2)
#		if ((separation_h[i]<match_tol) and our_h_halpha[i]!=99. and tendrup00_eqws[j]!=-9999.):
#			tendrup00_matches[i] = tendrup00_eqws[j]
#print tendrup00_matches

# let's just look at color/EqW for our MDM data
#good = np.where(our_h_halpha!=99.)[0]
#fig = plt.figure()
#ax = fig.add_subplot(1,1,1)
#plt.errorbar(our_h_color[good], our_h_halpha[good], our_h_halpha_err[good], capsize=0, ls='none', color='orange', marker='D')
#plt.axhline(y=0)
#plt.xlabel('($r^{\prime}-K$)')
#plt.ylabel('Measured EqW ($\AA$)')
#plt.xlim(1,5.25)
#plt.ylim(2.75,-10.)
#plt.title('Hyades')
##plt.show()

# look at the same for Stauffer 1991/1994/1997 (& Tendrup 2000)
#good = np.where(stauffer91_matches!=9999.)[0]
#good = np.where(stauffer97_matches!=9999.)[0]
##good = np.where(tendrup00_matches!=9999.)[0]
##print stauffer91_matches[good]
#fig = plt.figure()
#ax = fig.add_subplot(1,1,1)
#plt.plot(our_h_color[good], stauffer91_matches[good], ls='none', color='orange', marker='D')
#plt.plot(our_h_color[good], stauffer97_matches[good], ls='none', color='orange', marker='d')
#plt.plot(our_h_color[match], stauffer94_matches[match]*-1., ls='none', color='orange', marker='p')
##plt.plot(our_h_color[good], tendrup00_matches[good], ls='none', color='orange', marker='p')
#plt.axhline(y=0)
#plt.xlabel('($r^{\prime}-K$)')
#plt.ylabel('Stauffer et al. 1991, 1994, 1997 EqW ($\AA$)')
#plt.ylabel('Stauffer et al. 1997 EqW ($\AA$)')
##plt.ylabel('Stauffer et al. 1994 EqW ($\AA$)')
##plt.ylabel('Tendrup et al. 2000 EqW ($\AA$)')
#plt.xlim(1.5,5.5)
#plt.ylim(1.75,-6.0)
#plt.title('Hyades')
#plt.show()

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

# now want to look at median difference in EqW as a function of color
#good = np.where(our_h_halpha!=-99.)[0]
#our_good_halpha = our_h_halpha[good]
#clean = np.where(stauffer91_matches!=-9999.)
#clean_stauffer91 = stauffer91_matches[clean]
#diff_1 = (

#diff_h_err = np.sqrt(our_h_halpha_err[good]**2 + 0.026**2)
#diff_h_ha_1 = (our_h_halpha[good]-stauffer91_matches[good])
#diff_h_ha_2 = (our_h_halpha[good]-stauffer97_matches[good])
#diff_h_ha_3 = (our_h_halpha[match]-stauffer94_matches[match])
#clean = np.where(diff_h_ha_1!=10098.)[0]
#diff_1 = diff_h_ha_1[clean]
#color_1 = our_h_color[clean]
#err_1 = diff_h_err[clean]
#clean2 = np.where(diff_h_ha_2!=10098.)[0]
#diff_2 = diff_h_ha_1[clean2]
#color_2 = our_h_color[clean2]
#err_2 = diff_h_err[clean2]

#print diff_1

#for i in range(len(color_bins)-1):
#	in_this_bin = np.where((color_1>color_bins[i]) & 
#			(color_1<=color_bins[i+1]))[0]
#	binned_ha[i] = np.median(diff_1[in_this_bin]) #,diff_h_ha_2[in_this_bin],diff_h_ha_3[in_this_bin])
#	binned_errors[i] = np.median(err_1[in_this_bin])

#print binned_ha, binned_errors

#fig = plt.figure()
#ax = fig.add_subplot(1,1,1)
#lines = {'linestyle': 'None'}
#plt.plot(color_1, diff_1, color='orange', marker='D')
#plt.plot(our_h_color, diff_h_ha_2, color='orange', marker='d')
#plt.plot(our_h_color[match], diff_h_ha_3, color='orange', marker='p')
#plt.axhline(y=0, color='green')
#plt.errorbar(color_bins[:-1]+0.25, binned_ha, binned_errors, capsize = 0, ls='none', color ='red', marker = '*', ms=15)
#plt.xlabel('($r^{\prime}-K$)')
#plt.ylabel('Our EqW - Stauffer et al. 1991, 1994, 1997 EqW ($\AA$)')
#plt.xlim(-8.5,0.5)
#plt.xlim(1,6.0)
#plt.ylim(-3.75,1.5)
#plt.title('Hyades')
#plt.show()
#plt.savefig('eqw_praesepe_diff.eps', bbox_inches='tight')
