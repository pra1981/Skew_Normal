#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pylab import *
from scipy.special import erf
from scipy.stats import norm
from scipy.optimize import leastsq

import matplotlib.cm as cm
import matplotlib.colors as colors

sys.path.append('/Users/xavierdumusque/uni/phd/program/functions/')
sys.path.append('/Users/xavierdumusque/uni/phd/program/convective_blueshift/')

import glob,os,sys,string
import pyrdb
import functions
import pyfits
import pandas
import planet
import fit2

########################################################
#para
########################################################

stars = ['LRa01_E2_0165','HD192310','HD215152','HD128621']

path_dbase = '/Users/xavierdumusque/uni/phd/data/dbase/stars/'
path_data_Umberto = '../data/'
path_data_Pedro = 'pedrofigueira-line-profile-indicators/results/'

# From the Pepe paper, we only have lambda, the mean longitude
# Time of passage at periastron = epoch - (lambda - omega)/360*P = epoch - M/360*P = 2455151.02574596-(340.8-173)/360*74.72, where M is the mean anomalie
#planets = {'LRa01_E2_0165':[[3.42,0.85359165,4398.21,0.12,160*pi/180],[6.01,3.70,55953.3,0.12,160*pi/180]],'HD192310':[[3,74.72,55116.20,0.13,173*pi/180.],[2.27,525.8,55312,0.32,110*pi/180.]],'HD215152':[],'HD128621':[]}
planets = {'LRa01_E2_0165':[],'HD192310':[],'HD215152':[],'HD128621':[]}

########################################################
#read FITS files in 2010
########################################################

for star in stars:

    data_Umberto = pandas.read_csv(path_data_Umberto+star+'_Results_Umberto.dat',sep=' ',skiprows=[0],names=['index','jdb','vrad','bis_span','fwhm'])
    data_DRS = pandas.read_csv(path_dbase+star+'_harps_extract.rdb',sep='\t',skiprows=[1])
    file_Pedro = glob.glob(path_data_Pedro+star+'/WrapUp_*.txt')[0]
    data_Pedro = pandas.read_csv(file_Pedro,sep='\t\t',skiprows=[0,1],names=['file_root','jdb','vrad','bis_span','bis-','bis+','biGauss','Vasy','Vspan','fwhm'],engine='python')

    if star != 'HD128621':
        data_Umberto['jdb'] -= 2400000
    data_Pedro['jdb'] -= 2400000
    data_DRS['jdb'] = data_DRS['jdb'].round(6)
    data_Umberto['jdb'] = data_Umberto['jdb'].round(6)
    data_Pedro['jdb'] = data_Pedro['jdb'].round(6)
    if star != 'HD128621':
        data_Umberto[['vrad','bis_span','fwhm']] *= 1000.
    else:
        data_Umberto[['bis_span','fwhm']] *= 1000.
    data_DRS[['vrad','bis_span','fwhm']] *= 1000.
    data_Pedro[['vrad','bis_span','fwhm']] *= 1000.
    
    if star == 'HD128621':
        data_DRS['vrad'] -= (-22700.1747 - 0.5307 * (data_DRS['jdb'] - 55279.109840075726) - 1.83e-5 * (data_DRS['jdb'] - 55279.109840075726)**2)
        data_Pedro['vrad'] -= (-22700.1747 - 0.5307 * (data_Pedro['jdb'] - 55279.109840075726) - 1.83e-5 * (data_Pedro['jdb'] - 55279.109840075726)**2)
    
    data = data_DRS.merge(data_Umberto,on='jdb',how='inner')
    data = data.merge(data_Pedro,on='jdb',how='inner')
    data = data.loc[data['jdb'] < 57174.5] #remove new fibers

    data_DRS = data_DRS.loc[data_DRS['jdb'].isin(data['jdb'])]
    data_Umberto = data_Umberto.loc[data_Umberto['jdb'].isin(data['jdb'])]
    data_Pedro = data_Pedro.loc[data_Pedro['jdb'].isin(data['jdb'])]

    data_Umberto[['vrad','bis_span','fwhm']] -= data_Umberto[['vrad','bis_span','fwhm']].median()
    data_DRS[['vrad','bis_span','fwhm']] -= data_DRS[['vrad','bis_span','fwhm']].median()
    data_Pedro[['vrad','bis_span','fwhm']] -= data_Pedro[['vrad','bis_span','fwhm']].median()

    if star=='HD128621':
        data_Umberto = data_Umberto.loc[abs(data_Umberto['vrad']) < 20]
    elif star == 'HD215152':
        data_Umberto = data_Umberto.loc[abs(data_Umberto['vrad']) < 8]
    elif star=='HD192310':
        data_Umberto = data_Umberto.loc[abs(data_Umberto['vrad']) < 20]
        data_DRS = data_DRS.loc[abs(data_DRS['vrad']) < 20]

    data = data_DRS.merge(data_Umberto,on='jdb',how='inner')
    data = data.merge(data_Pedro,on='jdb',how='inner')

    data_DRS = data_DRS.loc[data_DRS['jdb'].isin(data['jdb'])]
    data_Umberto = data_Umberto.loc[data_Umberto['jdb'].isin(data['jdb'])]
    data_Pedro = data_Pedro.loc[data_Pedro['jdb'].isin(data['jdb'])]

    vrad_planet = zeros(len(data))
    for i,pl in enumerate(planets[star]):
        vrad_planet += planet.get_vrad_planet(data_DRS['jdb'],0,pl[0],pl[1],pl[2],time='periastron',e=pl[3],omega=pl[4],output='0',step_jdb_continuous=1.)

    data_DRS['vrad'] -= vrad_planet
    data_Umberto['vrad'] -= vrad_planet
    data_Pedro['vrad'] -= vrad_planet

    rcParams['figure.figsize'] = [15,8]

    figure()
    subplot(311)
    title(star)
    errorbar(data_DRS['jdb'],data_DRS['vrad'],data_DRS['svrad'],marker='o',ls='',zorder=0)
    scatter(data_DRS['jdb'],data_DRS['vrad'],c=data_DRS['sn50'],zorder=1)
    ylabel('DRS RV [m/s]')
    cb=colorbar()
    cb.set_label('SNR order 50')
    ax=gca()
    subplot(312,sharex=ax)
    errorbar(data_Umberto['jdb'],data_Umberto['vrad'],marker='o',ls='',zorder=0)
    scatter(data_Umberto['jdb'],data_Umberto['vrad'],c=data_DRS['sn50'],zorder=1)
    ylabel('SN RV [m/s]')
    xlabel('JD - 2400000 [d]')
    cb=colorbar()
    cb.set_label('SNR order 50')
    subplot(313)
    perio = functions.compute_periodogram_new(data_DRS['jdb'],data_DRS['vrad'],data_DRS['svrad'],ofac=4,min_P=1.1)
    semilogx(1./perio[0],perio[1])
    ylabel('Normalized Power')
    xlabel('Period [d]')

#    savefig('../figures/RV_vs_time_%s.pdf' % star)

    figure()
    para = ['vrad','bis_span']
    a,b = data_Umberto[para[0]],data_Umberto[para[1]]
    c,d = data_DRS[para[0]],data_DRS[para[1]]
    e,f = data_Pedro[para[0]],data_Pedro['bis-']
    g,h = data_Pedro[para[0]],data_Pedro['bis+']
    i,j = data_Pedro[para[0]],data_Pedro['biGauss']
    k,l = data_Pedro[para[0]],data_Pedro['Vasy']
    m,n = data_Pedro[para[0]],data_Pedro['Vspan']
    o,p = data_Pedro[para[0]],data_Pedro[para[1]]
    subplot(331)
    plot(a,b,'o',label='SNR=%.2f' % corrcoef(a,b)[0][1])
    ylabel('SN BIS [m/s]')
    xlabel('SN RV [m/s]')
    legend()
    subplot(332)
    title(star)
    plot(c,d,'o',label='SNR=%.2f' % corrcoef(c,d)[0][1])
    ylabel('DRS BIS [m/s]')
    xlabel('RV [m/s]')
    legend()
    subplot(333)
    plot(e,f,'o',label='SNR=%.2f' % corrcoef(e,f)[0][1])
    ylabel('FIGUEIRA BIS- [m/s]')
    xlabel('RV [m/s]')
    legend()
    subplot(334)
    ylabel('FIGUEIRA BIS+ [m/s]')
    xlabel('RV [m/s]')
    plot(g,h,'o',label='SNR=%.2f' % corrcoef(g,h)[0][1])
    legend()
    subplot(335)
    ylabel('FIGUEIRA biGauss [m/s]')
    plot(i,j,'o',label='SNR=%.2f' % corrcoef(i,j)[0][1])
    xlabel('RV [m/s]')
    legend()
    subplot(336)
    ylabel('FIGUEIRA Vasy [m/s]')
    plot(k,l,'o',label='SNR=%.2f' % corrcoef(k,l)[0][1])
    xlabel('RV [m/s]')
    legend()
    subplot(337)
    ylabel('FIGUEIRA Vspan [m/s]')
    plot(m,n,'o',label='SNR=%.2f' % corrcoef(m,n)[0][1])
    xlabel('RV [m/s]')
    legend()

    para = ['vrad','fwhm']
    a,b = data_Umberto[para[0]],data_Umberto[para[1]]
    c,d = data_DRS[para[0]],data_DRS[para[1]]
    e,f = data_Pedro[para[0]],data_Pedro[para[1]]
    subplot(338)
    plot(a,b,'o',label='SNR=%.2f' % corrcoef(a,b)[0][1],color='r')
    ylabel('SN FWHM [m/s]')
    xlabel('SN RV [m/s]')
    legend()
    subplot(339)
    plot(c,d,'o',label='SNR=%.2f' % corrcoef(c,d)[0][1],color='r')
    ylabel('DRS FWHM [m/s]')
    xlabel('RV [m/s]')
    legend()

    subplots_adjust(hspace=0.4,wspace=0.4)

    savefig('../figures/Comparison_para_%s.pdf' % star)


    Matrice = dstack([data_Umberto['bis_span'],data_Umberto['fwhm']])
    para1,sig_para1,mod1,chisq1 = fit2.linlsq(Matrice[0],array(data_Umberto['vrad']),array(data_DRS['svrad']))
    Umberto_activity_corr = sum([Matrice[0][:,i] * para1[i] for i in arange(len(para1))],0)

    Matrice = dstack([data_DRS['bis_span'],data_DRS['fwhm']])
    para1,sig_para1,mod1,chisq1 = fit2.linlsq(Matrice[0],array(data_DRS['vrad']),array(data_DRS['svrad']))
    DRS_activity_corr = sum([Matrice[0][:,i] * para1[i] for i in arange(len(para1))],0)

    figure()
    subplot(221)
    title(star +', SN RV : std=%.2f m/s' % std(data_Umberto['vrad']))
    plot(data_Umberto['jdb'],data_Umberto['vrad'],'o')
    ylabel('RV [m/s]')
    subplot(223)
    title('SN RV activity corr (SN RV - x*SN BIS - y*SN FWHM): std=%.2f m/s' % std(data_Umberto['vrad']-Umberto_activity_corr))
    plot(data_Umberto['jdb'],data_Umberto['vrad']-Umberto_activity_corr,'o')
    ylabel('RV [m/s]')
    xlabel('JD - 2400000 [d]')

    subplot(222)
    title('RV : std=%.2f m/s' % std(data_DRS['vrad']))
    plot(data_DRS['jdb'],data_DRS['vrad'],'o')
    subplot(224)
    title('RV activity corr (RV - x*BIS - y*FWHM): std=%.2f m/s' % std(data_DRS['vrad']-DRS_activity_corr))
    plot(data_DRS['jdb'],data_DRS['vrad']-DRS_activity_corr,'o')
    xlabel('JD - 2400000 [d]')

    savefig('../figures/Correction activity_%s.pdf' % star)


    figure()
    subplot(121)
    hist(data_Umberto['bis_span'],bins=100,label='SN')
    hist(data_DRS['bis_span'],bins=100,alpha=0.5,label='Normal')
    xlabel('bis span')
    legend()
    subplot(122)
    hist(data_Umberto['fwhm'],bins=100,label='SN')
    hist(data_DRS['fwhm'],bins=100,alpha=0.5,label='Normal')
    xlabel('fwhm')
    legend()

#   # Tried to remove the difference in velocity between Normal and Skew Normal, but it does not work
#    Matrice = dstack([array(data_Umberto['vrad']-data_DRS['vrad'])])
#    para1,sig_para1,mod1,chisq1 = fit2.linlsq(Matrice[0],array(data_DRS['vrad']),array(data_DRS['svrad']))
#    new_activity_corr = sum([Matrice[0][:,i] * para1[i] for i in arange(len(para1))],0)
#
#    figure()
#    subplot(211)
#    title('RV : std=%.2f m/s' % std(data_DRS['vrad']))
#    plot(data_DRS['jdb'],data_DRS['vrad'],'o')
#    ylabel('RV [m/s]')
#    subplot(212)
#    title('RV activity corr (SN RV - x*(SN_RV - RV)): std=%.2f m/s' % std(data_DRS['vrad']-new_activity_corr))
#    plot(data_DRS['jdb'],data_DRS['vrad']-new_activity_corr,'o')
#    ylabel('RV [m/s]')
#    xlabel('JD - 2400000 [d]')






