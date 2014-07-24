
from scipy.interpolate import interp1d
import numpy as np
import praesepe_comp as pc
import get_data, pickle, read_spec, bol_corr, get_sdss, triangle
import asciitable as at
import plot_grid as pg
from ha_cont import ha_cont
from emissionline import emissionline

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import Normalize
from matplotlib.patches import Ellipse
import matplotlib.cm as cm

kh = at.read('/home/stephanie/ptf/models/kraushillenbrand5.dat')
ksub = np.array([-23,-19,-15,-12,-11,-9,-7,-6])#arange(-15,-1,1)
kh_rpmK = (kh['Mr'] - 0.035*(kh['Mr']-kh['Mi']) - 0.007 - kh['MK'])[ksub]
kh_spt = kh['SpT'][ksub]
klen = len(kh_spt)

def hist2d(x, y, xbins=None, ybins=None, *args, **kwargs):
    """
    Plot a 2-D histogram of samples.
    Written by Dan Foreman-Mackey et al. (triangle.py)
    modified by Stephanie Douglas

    """
    ax = kwargs.pop("ax", plt.gca())

    extent = kwargs.pop("extent", [[x.min(), x.max()], [y.min(), y.max()]])
    bins = kwargs.pop("bins", 50)
    color = kwargs.pop("color", "k")
    linewidths = kwargs.pop("linewidths", None)
    plot_datapoints = kwargs.get("plot_datapoints", True)
    plot_contours = kwargs.get("plot_contours", True)

    cmap = cm.get_cmap("gray")
    cmap._init()
    cmap._lut[:-3, :-1] = 0.
    cmap._lut[:-3, -1] = np.linspace(1, 0, cmap.N)

    if xbins==None:
        X = np.linspace(extent[0][0], extent[0][1], bins + 1)
    else:
        X = xbins
    if ybins==None:
        Y = np.linspace(extent[1][0], extent[1][1], bins + 1)
    else:
        Y = ybins

    try:
        H, X, Y = np.histogram2d(x.flatten(), y.flatten(), bins=(X, Y))
    except ValueError:
        raise ValueError("It looks like at least one of your sample columns "
                         "have no dynamic range. You could try using the "
                         "`extent` argument.")

    V = 1.0 - np.exp(-0.5 * np.arange(0.5, 2.1, 0.5) ** 2)
    Hflat = H.flatten()
    inds = np.argsort(Hflat)[::-1]
    Hflat = Hflat[inds]
    sm = np.cumsum(Hflat)
    sm /= sm[-1]

    for i, v0 in enumerate(V):
        try:
            V[i] = Hflat[sm <= v0][-1]
        except:
            V[i] = Hflat[0]

    X1, Y1 = 0.5 * (X[1:] + X[:-1]), 0.5 * (Y[1:] + Y[:-1])
    X, Y = X[:-1], Y[:-1]

    if plot_datapoints:
        ax.plot(x, y, "o", color=color, ms=1.5, zorder=-1, alpha=0.1,
                rasterized=True)
#        if plot_contours:
#            ax.contourf(X1, Y1, H.T, [V[-1], H.max()],
#                        cmap=LinearSegmentedColormap.from_list("cmap",
#                                                               ([1] * 3,
#                                                                [1] * 3),
#                        N=2), antialiased=False)
#
#    if plot_contours:
#        ax.pcolor(X, Y, H.max() - H.T, cmap=cmap)
#        ax.contour(X1, Y1, H.T, V, colors=color, linewidths=linewidths)

    data = np.vstack([x, y])
    mu = np.mean(data, axis=1)
    cov = np.cov(data)
    if kwargs.pop("plot_ellipse", False):
        error_ellipse(mu, cov, ax=ax, edgecolor="r", ls="dashed")

    ax.set_xlim(extent[0])
    ax.set_ylim(extent[1])


pdat,pobs,pobs_nr,pobs_r = get_data.get_data('P')
hdat,hobs,hobs_nr,hobs_r = get_data.get_data('H')


colorbins = np.logspace(np.log10(0.5),np.log10(8),16)
colorcenter  = (colorbins[1:]+colorbins[:-1])/2.0
num_bins = len(colorcenter)
coloredges = np.zeros(num_bins*2).reshape(2,-1)
coloredges[0] = colorbins[1:]-colorcenter
coloredges[1] = colorcenter-colorbins[:-1]

hpmem = hdat.field('ROESER_PMEM')
ppmem = pdat.field('ADAMPMEM')
pmem_threshold = 70.0

p_avg = np.zeros(len(colorbins)-1)
h_avg = np.zeros(len(colorbins)-1)
m_avg = np.zeros(len(colorbins)-1)


p_std = np.zeros(len(colorbins)-1)
h_std = np.zeros(len(colorbins)-1)
p_std2 = np.zeros(len(colorbins)-1)
h_std2 = np.zeros(len(colorbins)-1)
m_std = np.zeros(len(colorbins)-1)

peqw = pdat.field('AVG_EQW')*-1.0
heqw = hdat.field('AVG_EQW')*-1.0
peqw_err = pdat.field('AVG_EQW_ERR')
heqw_err = hdat.field('AVG_EQW_ERR')

mchi, meqw, mlhalbol, mrpmK, mspt = get_sdss.get_halpha()

plt.figure(figsize=(9,8))
for i in range(len(colorbins)-1):
    pbin = ((pdat.field('RPRIME_K')>colorbins[i]) &
        (pdat.field('RPRIME_K')<=colorbins[i+1]) 
        & ((ppmem>=pmem_threshold) | (ppmem<0))
        #& (pbinary==False) 
        & (peqw<90))
    hbin = ((hdat.field('RPRIME_K')>colorbins[i]) &
        (hdat.field('RPRIME_K')<=colorbins[i+1])
        #& (hbinary==False) 
        & ((hpmem>=pmem_threshold) | (hpmem<0))
        & (heqw<90))
    mbin = ((mrpmK>colorbins[i]) & (mrpmK<=colorbins[i+1]))
    #print colorcenter[i]
    #print peqw[pbin]
    #print heqw[hbin]

    print colorcenter[i],len(np.where(pbin==True)[0]),
    print len(np.where(hbin==True)[0])

#    if (i<(len(colorbins)-2)):
#        plt.text(colorcenter[i],-10,str(len(np.where(pbin==True)[0])),
#            color='b')
#        plt.text(colorcenter[i],-9.5,str(len(np.where(hbin==True)[0])),
#            color='OrangeRed')

    p_avg[i] = np.average(peqw[pbin])
    p_std2[i] = np.sqrt(np.sum(peqw_err[pbin]*peqw_err[pbin]))/len(pbin)
    p_std[i] = np.std(peqw[pbin])
    h_avg[i] = np.average(heqw[hbin])
    h_std2[i] = np.sqrt(np.sum(heqw_err[hbin]*heqw_err[hbin]))/len(hbin)
    h_std[i] = np.std(heqw[hbin])

    m_avg[i] = np.average(meqw[mbin])
    m_std[i] = np.std(meqw[mbin])
    #print colorbins[i:i+2]
    #print heqw_err[hbin]



eqw_bins = np.arange(-15,10,0.5)
col_bins = colorbins
zeqw = np.histogram2d(meqw,mrpmK,(eqw_bins,col_bins))
logbins = np.arange(25,650,100)
norm = Normalize(-50,max(zeqw[0].flatten()))

inact_sum = 0
non_hist_sum = 0
for i in range(len(eqw_bins)-1):
    for j in range(len(col_bins)-1):
        if (zeqw[0][i,j]>24) and (eqw_bins[i]>=-3):
            inact_sum += zeqw[0][i,j]
        else:
            non_hist_sum += zeqw[0][i,j]
#            if zeqw[0][i,j]>0:
#                print "not histogram - ", zeqw[0][i,j]
print "{} stars in the 'inactive' area of the histogram".format(inact_sum)
print "(Does not include that one active patch at late types)"
print "{} stars that are individual points, and that one 25 star patch".format(
    non_hist_sum)
print "and {} SDSS stars that have EqW <-15, so aren't included".format(
    len(np.where(meqw<-15)[0]))
print "{} total stars by adding, {} in the sample".format(
    inact_sum+non_hist_sum+len(np.where(meqw<-15)[0]),len(meqw))


z2 = zeqw[0].flatten()
lloc = np.where(z2<25)[0]
z2[lloc] = -50.
for l in logbins:
    #print z2,l,l-100
    lloc = np.where((z2>=l) & (z2<(l+100)))[0]
    z2[lloc] = l
z3 = z2.reshape(zeqw[0].shape)



histogram = plt.pcolor(col_bins,eqw_bins,z3,cmap=cm.get_cmap("Greys"),norm=norm)

for i in range(len(eqw_bins)-1):
    for j in range(len(col_bins)-1):
        if z3[i,j]>24:
            bad = np.where((meqw>=eqw_bins[i]) & (meqw<eqw_bins[i+1]) &
               (mrpmK>=col_bins[j]) & (mrpmK<col_bins[j+1]))[0]
            meqw = np.delete(meqw,bad)
            mrpmK = np.delete(mrpmK,bad)

#Need to make these points show up now
plt.plot(mrpmK,meqw,'o',color='#C0C0C0', ms=1.5,mec='none')
#    , zorder=100, alpha=0.1, rasterized=True)

#eqw_bins = np.delete(eqw_bins,-1)+0.25
#col_bins = colorcenter
#contour_fill = plt.contourf(col_bins,eqw_bins,zeqw[0],
#     logbins,cmap=cm.get_cmap("Greys"),norm=norm)
#contour_lines = plt.contour(col_bins,eqw_bins,zeqw[0],
#     logbins,colors='Grey',label='SDSS Field')
#print logbins

hist2d(mrpmK,meqw,col_bins,bins=100)
triangle.hist2d(mrpmK,meqw,plot_contours=False)

plt.errorbar(colorcenter,p_avg,p_std,coloredges,fmt='o',color='b',
    mec='b',label='Praesepe',
    capsize=0)
plt.errorbar(colorcenter+0.01,h_avg,h_std,coloredges,fmt='D',mfc='OrangeRed',
    color='OrangeRed',mec='OrangeRed',label='Hyades',capsize=0)
plt.xlabel('(r\'-K)',fontsize='x-large')
plt.ylabel(r'Average $H\alpha$ EqW',fontsize='x-large')
plt.ylim(7,-11)
plt.xlim(0.25,6.8)
plt.legend(loc=4,numpoints=1)

#plt.clabel(contour_lines,fmt='%i',fontsize=10,manual=True,rightside_up=True)
ax = plt.gca()
ax.tick_params(which='both',labelsize='large',top=False)
texty = -11.25
for i in range(klen):
   ax.text(kh_rpmK[i],texty,kh_spt[i],fontsize='large')

plt.show()

plt.savefig('papereqws_compare.png')
plt.savefig('papereqws_compare.ps')

