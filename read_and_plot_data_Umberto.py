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

########################################################
#para
########################################################

stars = ['LRa01_E2_0165','HD192310','HD215152','HD128621','HD10700']

path_dbase = '/Users/xavierdumusque/uni/phd/data/dbase/stars/'
path_data_Umberto = '../data/skew_normal/'
path_data_Pedro = 'pedrofigueira-line-profile-indicators/results/'

########################################################
#read FITS files in 2010
########################################################

for star in stars:

    data_Umberto = pandas.read_csv(path_data_Umberto+star+'_Results_Umberto.dat',sep=' ',skiprows=[0],names=['index','jdb','vrad','bis_span','fwhm'])
    data_DRS = pandas.read_csv(path_dbase+star+'_harps_extract.rdb',sep='\t',skiprows=[1])
    file_Pedro = glob.glob(path_data_Pedro+star+'/WrapUp_*.txt')[0]
    data_Pedro = pandas.read_csv(file_Pedro,sep='\t\t',skiprows=[0,1],names=['file_root','jdb','vrad','bis_span','bis-','bis+','biGauss','Vasy','Vspan','fwhm'],engine='python')

    data_Umberto['jdb'] -= 2400000
    data_Pedro['jdb'] -= 2400000
    data_Umberto[['vrad','bis_span','fwhm']] *= 1000.
    data_DRS[['vrad','bis_span','fwhm']] *= 1000.
    data_Pedro[['vrad','bis_span','fwhm']] *= 1000.

    data_Umberto[['vrad','bis_span','fwhm']] -= data_Umberto[['vrad','bis_span','fwhm']].mean()
    data_DRS[['vrad','bis_span','fwhm']] -= data_DRS[['vrad','bis_span','fwhm']].mean()
    data_Pedro[['vrad','bis_span','fwhm']] -= data_Pedro[['vrad','bis_span','fwhm']].mean()
    
    figure()
    subplot(311)
    plot(data_DRS['jdb'],data_DRS['vrad'])
    ax=gca()
    subplot(312,sharex=ax)
    plot(data_Umberto['jdb'],data_Umberto['vrad'])
    subplot(313,sharex=ax)
    plot(data_Pedro['jdb'],data_Pedro['vrad'])

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
    title('SN BIS %s' % star)
    plot(a,b,'o',label='SN R=%.2f' % corrcoef(a,b)[0][1])
    legend()
    subplot(332)
    title('DRS BIS')
    plot(c,d,'o',label='SN R=%.2f' % corrcoef(c,d)[0][1])
    legend()
    subplot(333)
    title('FIGUEIRA BIS-')
    plot(e,f,'o',label='SN R=%.2f' % corrcoef(e,f)[0][1])
    legend()
    subplot(334)
    title('FIGUEIRA BIS+')
    plot(g,h,'o',label='SN R=%.2f' % corrcoef(g,h)[0][1])
    legend()
    subplot(335)
    title('FIGUEIRA biGauss')
    plot(i,j,'o',label='SN R=%.2f' % corrcoef(i,j)[0][1])
    legend()
    subplot(336)
    title('FIGUEIRA Vasy')
    plot(k,l,'o',label='SN R=%.2f' % corrcoef(k,l)[0][1])
    legend()
    subplot(337)
    title('FIGUEIRA Vspan')
    plot(m,n,'o',label='SN R=%.2f' % corrcoef(m,n)[0][1])
    legend()
#    subplot(338)
#    title('FIGUEIRA BIS')
#    plot(o,p,'o',label='SN R=%.2f' % corrcoef(o,p)[0][1])
#    legend()

    savefig('../figures/skew_normal/Comparison_BIS_%s.pdf' % star)

    figure()
    para = ['vrad','fwhm']
    a,b = data_Umberto[para[0]],data_Umberto[para[1]]
    c,d = data_DRS[para[0]],data_DRS[para[1]]
    e,f = data_Pedro[para[0]],data_Pedro[para[1]]
    subplot(131)
    title('SN FWHM %s' % star)
    plot(a,b,'o',label='SN R=%.2f' % corrcoef(a,b)[0][1])
    legend()
    subplot(132)
    title('DRS FWHM')
    plot(c,d,'o',label='SN R=%.2f' % corrcoef(c,d)[0][1])
    legend()
    subplot(133)
    title('FIGUEIRA FWHM')
    plot(e,f,'o',label='SN R=%.2f' % corrcoef(e,f)[0][1])
    legend()

    savefig('../skew_normal/figures/Comparison_FWHM_%s.pdf' % star)




