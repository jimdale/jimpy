#!/usr/bin/python2.7

import os
import sys

import numpy
from numpy import *

from scipy.stats import ks_2samp


import healpy
import string

import matplotlib
import matplotlib
matplotlib.use('Agg')

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
from pylab import *
from matplotlib import colors
from matplotlib import rc
rc('font',**{'family':'serif','serif':['Times New Roman'],'size':12})
rc('text',usetex=True)
rc('patch',antialiased=False)

from scipy.stats import distributions

# https://github.com/keflavich/plfit
import plfit

import moment_maps_rotate

def analyze(fname, xmin=-20., xmax=20., ymin=-20., ymax=20., zmin=-100.,
            zmax=100., limitsmean=[0.,2,-2.,6.], iline=100):
    #variables for screenshots - NB imcol expects units in CODE UNITS

    subplots_adjust(hspace=0.1)
    subplots_adjust(wspace=0.1)

    # conversion factor gcm^-2 --> Av
    extfactor=1./1.67e-24/1.e21

    # physical limits on field-of-view
    # [xyz]minmax

    # pixel dimensions of images
    #iline=100

    # dimensions of plot array

    nplotx=6
    nploty=4

    # size of multi-plots

    fsize=(18,27)

    # column-density colorbar limits

    logsigmax=0.3
    logrange=4
    logsigmin=logsigmax-logrange

    # velocity colorbar limits

    vmax=5.
    vmin=-5.

    # limits for mean pdf:
    #limitsmean=[0.,2,-2.,6.]

    # lower limit for cdf column densities
    #cmflowlimit=limitsmean[0]
    cmflowlimit=0.
    # set negative for no cmfhighlimit
    cmfhighlimit=-1.


    normsig=colors.Normalize(vmin=logsigmin,vmax=logsigmax)
    normvel=colors.Normalize(vmin=vmin,vmax=vmax)
    normsigav=colors.Normalize(vmin=log10((10.**logsigmin)*extfactor),vmax=log10((10.**logsigmax)*extfactor))

    # temperature thresholds (gas below coldtemp and above hottemp is rejected)

    coldtemp=0.
    hottemp=125.



    # masks for non-antiparallel Healpix vectors

    level0mask=[0,1,2,3,4,5]
    level1mask=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23]

    # set Healpix level

    nhplevel=1
    nhppix=12*4**nhplevel

    if nhplevel==0:
        mask=level0mask
    elif nhplevel==1:
        mask=level1mask
    else:
        print 'This script is too dumb to handle Healpix level',nhplevel
        raise IndexError

    print 'Healpix level: ',nhplevel
    print 'Healpix wedges: ',nhppix

    healpix_angles = zip(*healpy.pix2ang(nhplevel+1, np.arange(nhppix)))
    healpix_vectors = zip(*healpy.pix2vec(nhplevel+1, np.arange(nhppix)))
    angular_separations = np.array([[np.arccos(np.dot(x,y))
                                     for x in healpix_vectors[:nhppix/2]]
                                    for y in healpix_vectors[:nhppix/2]])

    rad_to_deg=360./(2.*pi)

    xlab='x(pc)'
    ylab='y(pc)'

    #plt.tight_layout()

    shr=0.73

    extent=(xmin,xmax,ymin,ymax)


    moment0maps=[]
    moment1maps=[]
    pdfs=[]
    data = None

    for ihppix in range(nhppix):
        if (ihppix) in mask:
            print 'Plotting: ',ihppix,' of ',nhppix,' images, nplotx,y=',nplotx,nploty
         
            angles=healpix_angles[ihppix]
         
            angle1=angles[1]*rad_to_deg
            angle2=0.0
            angle3=angles[0]*rad_to_deg
         
         
            print 'Angles: ',angle1,angle2,angle3
         
             #call wrapper for imcol
         
            fnamei=fname
         
            coldens, vels = moment_maps_rotate.imager(fnamei, xmin, xmax, ymin,
                                                      ymax, zmin, zmax, logsigmax,
                                                      logrange, logsigmin, angle1,
                                                      angle2, angle3, iline,
                                                      coldtemp, hottemp)
         
            moment0 = coldens
            moment1 = vels
            
         
            moment0maps.append(moment0)
            moment1maps.append(moment1)


    nmaps=len(moment0maps)

    ##########################################################################################
    #
    #   ZEROTH MOMENT MAPS
    #
    ##########################################################################################

    fig=plt.figure(1,figsize=fsize)
    clf()
    cla()

    set_cmap('hot')

    nplot=0

    for i in range(nmaps):
        nplot+=1
        plt.subplot(nplotx,nploty,nplot)

        moment0=moment0maps[nplot-1]

        imgplot=plt.imshow(moment0+np.log10(extfactor),extent=extent,origin='lower',interpolation='nearest')
        # imgplot.set_norm(normsig)
        imgplot.set_norm(normsigav)
                
        plt.contour(moment0*extfactor,levels=[10.**0.5,10.**1.],colors='white')
        
        plt.title('Wedge '+str(i))
        plt.xlim(xmin, xmax)
        plt.ylim(ymin, ymax)
        xlabel(xlab)
        ylabel(ylab)

    cax = plt.axes([0.95, 0.1, 0.03, 0.8])


    cb=plt.colorbar(imgplot, cax=cax)
    cb.set_label('Log(column density) (g cm$^{-2}$)')

    plt.savefig('moment0_angles_'+os.path.split(fname)[-1]+'.png',dpi=300,bbox_inches='tight')


    ##########################################################################################
    #
    #   FIRST MOMENT MAPS
    #
    ##########################################################################################

    fig=plt.figure(1,figsize=fsize)
    clf()
    cla()

    set_cmap('rainbow')

    nplot=0

    for i in range(nmaps):
            nplot+=1
            plt.subplot(nplotx,nploty,nplot)

            moment1=moment1maps[nplot-1]

            imgplot=plt.imshow(moment1,extent=extent,origin='lower',interpolation='nearest')
            imgplot.set_norm(normvel)

            plt.xlim(xmin, xmax)
            plt.ylim(ymin, ymax)
            xlabel(xlab)
            ylabel(ylab)

    cax = plt.axes([0.95, 0.1, 0.03, 0.8])

    cb=plt.colorbar(imgplot, cax=cax)
    cb.set_label('Line-of-sight velocity (km s$^{-1}$)')
    plt.savefig('moment1_angles_'+os.path.split(fname)[-1]+'.png',dpi=300,bbox_inches='tight')


    ##########################################################################################
    #
    #   MEAN PDF
    #
    ##########################################################################################

    fig=plt.figure(2,figsize=(8,8))
    clf()
    cla()

    for i in range(nmaps):
            pdf=ravel(moment0maps[i])
            pdf=[log10((10.**p)*extfactor) for p in pdf]
            pdfs.append(pdf)
    #       print pdf

    meanpdf=[]
    stdpdf=[]
    actualns=[]
    cmfxs=[]
    cmfys=[]
    nbins=50
    for coldenslist in pdfs:
            n,bins,patches=plt.hist(coldenslist,bins=nbins,range=(-3,3),log=True,histtype='step',fill=False)
            actualns.append(n)
    # crop cdf at lower surf dense limit
            cmf=sort(coldenslist)
            nlist=len(coldenslist)
            for i in range(nlist):
                    if cmf[i]>cmflowlimit:
                            cmftrim=cmf[i:]
                            ncmf=len(cmftrim)
                            break

                    
            if cmfhighlimit<0.:
                    ncmf=len(cmftrim)
                    fncmf=float(ncmf)
                    cmfxs.append(cmftrim)
                    cmfys.append((arange(ncmf))/fncmf)
            else:
                    for j in range(ncmf):
                            cmftrim=cmftrim[:j]
    #                       print i,j
                            ncmf=len(cmftrim)
                            fncmf=float(ncmf)
                            cmfxs.append(cmftrim)
                            cmfys.append((arange(ncmf))/fncmf)
                            break

    # perform KS tests of all against all

    fig=plt.figure(2,figsize=(8,8))
    clf()
    cla()

    for i in range(nmaps):
            for j in range(nmaps):
                    if i!=j:
                            ks=ks_2samp(cmfxs[i], cmfxs[j])
                            plt.plot(float(i+1),ks[0],'bx')
                            print i,j,ks
       
    xlabel('Healpix wedge')
    ylabel('KS statistics')

    plt.savefig('ks_stats_'+os.path.split(fname)[-1]+'.png',dpi=300,bbox_inches='tight')

    fig=plt.figure(2,figsize=(8,8))
    clf()
    cla()

    for i in range(nmaps):
            for j in range(nmaps):
                    if i!=j:
                            ks=ks_2samp(cmfxs[i], cmfxs[j])
                            plt.plot(float(i+1),ks[1],'bx')
       
    xlabel('Healpix wedge')
    ylabel('KS p-value')

    plt.savefig('ks_p-values_'+os.path.split(fname)[-1]+'.png',dpi=300,bbox_inches='tight')



    def ks_2samp_loc(data1,data2):
        data1 = np.sort(data1)
        data2 = np.sort(data2)
        n1 = data1.shape[0]
        n2 = data2.shape[0]
        data_all = np.concatenate([data1, data2])
        cdf1 = np.searchsorted(data1, data_all, side='right') / (1.0*n1)
        cdf2 = np.searchsorted(data2, data_all, side='right') / (1.0*n2)
        d_loc = np.argmax(np.absolute(cdf1 - cdf2))
        d = np.max(np.absolute(cdf1 - cdf2))
        # Note: d absolute not signed distance
        en = np.sqrt(n1 * n2 / float(n1 + n2))
        prob = distributions.kstwobign.sf((en + 0.12 + 0.11 / en) * d)

        return (d, prob, data_all[d_loc])


    fig=plt.figure(2,figsize=(8,8))
    clf()
    cla()

    grid = zeros((nmaps,nmaps))*np.nan
    grid_loc = zeros((nmaps,nmaps))*np.nan
    for i in range(nmaps):
            for j in range(nmaps):
                    if i<=j:
                            ks_D, ks_prob, ks_loc =ks_2samp_loc(cmfxs[i], cmfxs[j])
                            grid[i,j]=ks_prob
                            grid_loc[j,i]=ks_loc
       
    plt.imshow(grid,cmap='hot',interpolation='nearest',origin='lower')
    cb = colorbar()
    cb.set_label("KS $p$ value $P(D)$")
    plt.imshow(10**grid_loc,cmap='Blues',interpolation='nearest',origin='lower')
    cb = colorbar()
    cb.set_label("$A_V$")
    xlabel('Healpix wedge ID')
    ylabel('Healpix wedge ID')

    plt.savefig('ks_grid_'+os.path.split(fname)[-1]+'.png',dpi=300,bbox_inches='tight')

    plt.figure(2,figsize=(8,8))
    clf()
    cla()

    plt.hist(10**grid_loc.T[grid<0.05], color='r', histtype='step', label='Different ($p<0.05$)')
    plt.hist(10**grid_loc.T[grid>0.05], color='k', histtype='step', label='Same ($p>0.05$)')
    plt.legend(loc='best')
    xlabel("$A_V$ where the KS-test was computed")

    plt.savefig('ks_location_histograms_'+os.path.split(fname)[-1]+'.png',dpi=300,bbox_inches='tight')

    plt.figure(2,figsize=(8,8))
    clf()
    cla()

    plt.plot((angular_separations*180/np.pi).flat, grid.flat, 'b.')
    plt.plot((180-angular_separations*180/np.pi).flat, grid.flat, 'b.')
    plt.xlabel("Angular Separation ($^{\circ}$)")
    plt.ylabel("KS Probabilities")

    plt.savefig('ks_vs_angsep_'+os.path.split(fname)[-1]+'.png',dpi=300,bbox_inches='tight')

    ok = (grid!=1.0) & (angular_separations > 1e-4) & (np.isfinite(grid))
    same_ang_sep = np.concatenate([angular_separations[(grid>0.05) & ok],
                                   np.pi-angular_separations[(grid>0.05) & ok]])
    diff_ang_sep = np.concatenate([angular_separations[(grid<0.05) & ok],
                                   np.pi-angular_separations[(grid<0.05) & ok]])

    clf()
    cla()

    bins2 = np.linspace(0,90,10)
    plt.hist(diff_ang_sep*180/np.pi, color='r', histtype='step', label='Different ($p<0.05$)', bins=bins2)
    plt.hist(same_ang_sep*180/np.pi, color='k', histtype='step', label='Same ($p>0.05$)', bins=bins2)
    plt.hist((angular_separations[ok]*180/np.pi),
             color='b', histtype='step', label='Total # of images', bins=bins2)
    plt.legend(loc='upper left')
    plt.xlabel("Angular Separation ($^{\circ}$)")
    plt.xlim(0,90)

    plt.savefig('ks_vs_angsep_histograms_'+os.path.split(fname)[-1]+'.png',dpi=300,bbox_inches='tight')

    clf()
    cla()
    diff_counts,diff_edges = np.histogram(diff_ang_sep*180/np.pi, bins=bins2)
    same_counts,same_edges = np.histogram(same_ang_sep*180/np.pi, bins=bins2)
    total_counts_a = diff_counts+same_counts
    total_counts,total_edges = np.histogram(np.concatenate([(angular_separations[ok]*180/np.pi),
                                                            180-(angular_separations[ok]*180/np.pi)]),
                                            bins=bins2)

    xaxis = np.ravel(zip(diff_edges[:-1], diff_edges[1:]))
    yaxis = np.ravel(zip((diff_counts/total_counts), (diff_counts.astype('float')/total_counts)))
    plt.plot(xaxis, yaxis, drawstyle='steps',
            color='r', linewidth=2, alpha=0.7,
            label='Different ($p<0.05$)')
    yaxis = np.ravel(zip((same_counts/total_counts), (same_counts.astype('float')/total_counts)))
    plt.plot(xaxis, yaxis, drawstyle='steps',
             color='b', linewidth=2, alpha=0.7,
             label='Same ($p>0.05$)')

    plt.legend(loc='upper left')
    plt.xlabel("Angular Separation ($^{\circ}$)")
    plt.ylabel("Fraction of images")
    plt.xlim(0,90)

    plt.savefig('ks_vs_angsep_normalized_histograms_'+os.path.split(fname)[-1]+'.png',dpi=300,bbox_inches='tight')
    assert all(total_counts_a==total_counts)


    for j in range(nbins):
        binvalues=[]
        for i in range(nmaps):
            binvalues.append((actualns[i][j]))
        av=mean(binvalues)
        st=std(binvalues)
        if av<1e-2:av=1e-2
        meanpdf.append(av)
        stdpdf.append(st)

    high=[]
    low=[]
    for i in range(nbins):
        hi=meanpdf[i]+stdpdf[i]
        if hi<0:hi=1e-2
        high.append(hi)
        lo=meanpdf[i]-stdpdf[i]
        if lo<0:lo=1e-2
        low.append(lo)

    bns=[]
    for i in range(nbins):
            bns.append(0.5*(bins[i]+bins[i+1]))
            
    sh=sum(high)
    sm=sum(meanpdf)
    sl=sum(low)

    high=[log10(h) for h in high]
    meanpdf=[log10(m) for m in meanpdf]
    low=[log10(l) for l in low]

    plt.figure(2, figsize=(8,8))
    clf()
    cla()

    #plt.plot(bins,high,'r--')
    plt.plot(bns,meanpdf,'b-')
    #plt.plot(bins,meanpdfi,'r-')
    #plt.plot(bins,low,'r--')
    plt.fill_between(bns, low, high, edgecolor=None,facecolor='blue',alpha=0.5,interpolate=True)
    #plt.fill_between(bins, lowi, highi, edgecolor=None,facecolor='red',alpha=0.5,interpolate=True)
    xlab='log A$_{\\rm V}$'
    ylab='log P(A$_{\\rm V}$)'
    xlabel(xlab)
    ylabel(ylab)
    axis(limitsmean)
    plt.savefig('moment0_meanpdf'+os.path.split(fname)[-1]+'.png',dpi=300,bbox_inches='tight')

    ##########################################################################################
    #
    #   INDIVIDUAL CMFS
    #
    ##########################################################################################

    fig=plt.figure(1,figsize=fsize)
    clf()
    cla()

    nplot=0
    for i in range(nmaps):
    #       nplot+=1
    #        plt.subplot(nplotx,nploty,nplot)
            plt.plot(cmfxs[i],cmfys[i])
            xlab='log A$_{\\rm V}$'
            ylab='log P(A$_{\\rm V}$)'
            xlabel(xlab)
            ylabel(ylab)
    #       axis(limitsmean)

    plt.savefig('cmfs_angles'+os.path.split(fname)[-1]+'.png',dpi=300,bbox_inches='tight')

    ##########################################################################################
    #
    #   INDIVIDUAL PDFS COMPARED TO MEAN
    #
    ##########################################################################################

    fig=plt.figure(1,figsize=fsize)
    clf()
    cla()

    nplot=0
    for i in range(nmaps):
            nplot+=1
            plt.subplot(nplotx,nploty,nplot)
            plt.plot(bns,meanpdf,'b-')
            pdfi=actualns[i]
            pdfi=[log10(p) for p in pdfi]
            plt.plot(bns,pdfi,'k-')
            xlab='log A$_{\\rm V}$'
            ylab='log P(A$_{\\rm V}$)'
            xlabel(xlab)
            ylabel(ylab)
            axis(limitsmean)

    plt.savefig('pdfs_angles'+os.path.split(fname)[-1]+'.png',dpi=300,bbox_inches='tight')


    fig=plt.figure(1,figsize=fsize)
    clf()
    cla()

    nplot=0
    for i in range(nmaps):
            plt.plot(bns,meanpdf,'-')
            pdfi=actualns[i]
            pdfi=[log10(p) for p in pdfi]
            plt.plot(bns,pdfi,'k-')
            xlab='log A$_{\\rm V}$'
            ylab='log P(A$_{\\rm V}$)'
            xlabel(xlab)
            ylabel(ylab)
            axis(limitsmean)

    plt.savefig('all_pdfs_angles'+os.path.split(fname)[-1]+'.png',dpi=300,bbox_inches='tight')


    ks_pl = []
    ks_ln = []
    ks_pl_nofit = []
    ks_ln_nofit = []

    for row in cmfxs:
        try:
            pf = plfit.plfit(10**np.array(row), discrete=False, usefortran=True)
            pf.lognormal()
            ks_pl.append(pf._ks_prob)
            ks_ln.append(pf.lognormal_ksP)
        except AssertionError:
            print("There was an error. Skipping a value")

        try:
            pf = plfit.plfit(10**np.array(row), xmin=10**0.5, discrete=False, usefortran=True)
            pf.lognormal()
            ks_pl_nofit.append(pf._ks_prob)
            ks_ln_nofit.append(pf.lognormal_ksP)
        except AssertionError:
            print("There was an error. Skipping a value")

    plt.figure(2, figsize=(8,8))
    plt.clf()
    try:
        plt.plot(ks_pl, ks_ln, 'bo')
        plt.plot(ks_pl_nofit, ks_ln_nofit, 'ro')
    except:
        print("Failed to plot plfits")

    # plot cutoffs
    plt.plot([0.05,0.05],[0,1],'k--')
    plt.plot([0,1],[0.05,0.05],'k--')    
    plt.xlabel("Powerlaw KS Probability")
    plt.ylabel("Lognormal KS Probability")

    clf()
    cla()
    try:
        plt.plot(ks_pl_nofit, ks_pl, 'o')
        plt.plot([0,1],[0,1],'k--')
        plt.title("Everything should be above the line")
    except:
        print("Failed to plot plfits")


    print 'Done!'

if __name__ == "__main__":

    for ii in range(3):
        plt.close(ii)

    params = {'/Users/adam/work/jimsims/simulation_data/ext_feedback/BTOP100': {'xmin':-60., 'xmax':60., 'ymin':-60., 'ymax':60., 'limitsmean':[0.,3.,-2,2.5]},
              '/Users/adam/work/jimsims/simulation_data/ext_feedback/BTOP150': {'xmin':-60., 'xmax':60., 'ymin':-60., 'ymax':60., 'limitsmean':[0.,3.,-2,2.5]},
              '/Users/adam/work/jimsims/simulation_data/ext_feedback/BTOP200': {'xmin':-60., 'xmax':60., 'ymin':-60., 'ymax':60., 'limitsmean':[0.,3.,-2,2.5]},
              '/Users/adam/work/jimsims/simulation_data/ext_feedback/BTOP350': {'xmin':-60., 'xmax':60., 'ymin':-60., 'ymax':60., 'limitsmean':[0.,3.,-2,2.5]},
              '/Users/adam/work/jimsims/simulation_data/run_i/RUNI180':  {'xmin':-20., 'xmax':20., 'ymin':-20., 'ymax':20., 'limitsmean':[0.,2,-2.,3.]},}

    for k,v in params.items():
        analyze(k, **v)
