import triangle, cPickle
import numpy as np
import matplotlib.pyplot as plt


def quantile(x,quantiles):
    xsorted = sorted(x)
    qvalues = [xsorted[int(q * len(xsorted))] for q in quantiles]
    return zip(quantiles,qvalues)


#infile = open('/home/stephanie/ptf/fit_rossby/fit_rossby.pkl','rb')
infile = open('fit_rossby.pkl','rb')
samples = cPickle.load(infile)
infile.close()
sl_mcmc = quantile(samples[:,0],[.16,.5,.84]) # saturation level
to_mcmc = quantile(samples[:,1],[.16,.5,.84]) # turnover
be_mcmc = quantile(samples[:,2],[.16,.5,.84]) # beta

ndim=3

triangle.corner(samples.reshape((-1,ndim)))#,quantiles=[0.16,0.50,0.84])

# big labels
blabel = r'$\beta$'
rolabel = r'$Ro_{sat}$'
llabel = r'$(L_{H\alpha}/L_{bol})_{sat}\ (\times10^{-4})$'

ax4 = plt.subplot(334)
ax4.set_ylabel(rolabel,fontsize='x-large')

ax7 = plt.subplot(337)
xticks = ax7.get_xticks()
new_labels = []
for xt in xticks:
    new_labels =np.append(new_labels,str(xt*10000))
ax7.set_xticklabels(new_labels)
ax7.set_xlabel(llabel,fontsize='x-large')
ax7.set_ylabel(blabel,fontsize='x-large')

ax8 = plt.subplot(338)
ax8.set_xlabel(rolabel,fontsize='x-large')

ax9 = plt.subplot(339)
ax9.set_xlabel(blabel,fontsize='x-large')

plt.savefig('paper_corner.png',dpi=300,bbox_inches='tight')


# Add quantiles/lit values
ax1 = plt.subplot(331)
yl = ax1.get_ylim()
ax1.plot((sl_mcmc[1][1],sl_mcmc[1][1]),yl,'-',color='grey')
ax1.plot((sl_mcmc[0][1],sl_mcmc[0][1]),yl,'-',color='grey')
ax1.plot((sl_mcmc[2][1],sl_mcmc[2][1]),yl,'-',color='grey')


ax5 = plt.subplot(335)
yl = ax5.get_ylim()
ax5.plot((to_mcmc[1][1],to_mcmc[1][1]),yl,'-',color='grey')
ax5.plot((to_mcmc[0][1],to_mcmc[0][1]),yl,'-',color='grey')
ax5.plot((to_mcmc[2][1],to_mcmc[2][1]),yl,'-',color='grey')


yl = ax9.get_ylim()
ax9.plot((be_mcmc[1][1],be_mcmc[1][1]),yl,'-',color='grey')
ax9.plot((be_mcmc[0][1],be_mcmc[0][1]),yl,'-',color='grey')
ax9.plot((be_mcmc[2][1],be_mcmc[2][1]),yl,'-',color='grey')
#ax9.plot((-1,-1),yl,'k--') # Jackson & Jeffries

plt.savefig('paper_corner_quantiles.png',dpi=300,bbox_inches='tight')

