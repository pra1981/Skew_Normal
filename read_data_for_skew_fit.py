#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pylab import *
from scipy.special import erf
from scipy.stats import norm
from scipy.optimize import leastsq

import matplotlib.cm as cm
import matplotlib.colors as colors

sys.path.append('/Users/xav/uni/phd/program/functions/')
sys.path.append('/Users/xav/uni/phd/program/convective_blueshift/')
sys.path.append('/home/spectro/dumusque/phd/program/functions/')

import glob,os,sys,string
import pyrdb
import functions
import pyfits
import pandas as pd

########################################################
#para
########################################################

#dbase
dbase_folder = '/dace/data/observations/harps/DRS-3.5/dbase/'

#reduced files
reduced_folder = '/dace/data/observations/harps/DRS-3.5/reduced/'

star = raw_input('For which star do you want to extract the CCFs and other infos: ')

extract_file = dbase_folder+star+'_harps_extract.rdb'

########################################################
#read FITS files in 2010
########################################################

#Read good meas
data = pd.read_csv(extract_file,sep='\t',skiprows=[1])
len_data = len(data)


time,CCF,vrad,sn50,rv,bis,fwhm,contrast = [],[],[],[],[],[],[],[]
for i,(night,fileroot,mask) in enumerate(zip(data['night'],data['file_root'],data['mask'])):

    if i%100 == 0:
        print i, ' over ',len_data
    filename = reduced_folder+night+'/'+fileroot+'_ccf_'+mask+'_A.fits'
    header = pyfits.getheader(filename)
    time_tmp = header['ESO DRS BJD']

    time = append(time,time_tmp)
    CCF.append((pyfits.getdata(filename)[-1]).astype('float64'))
    vrad.append(arange(header['CRVAL1'],header['CRVAL1']+len(CCF[-1])*header['CDELT1'],header['CDELT1']))
    sn50 = append(sn50,header['ESO DRS SPE EXT SN50'])

    header_bis = pyfits.getheader(reduced_folder+night+'/'+fileroot+'_bis_'+mask+'_A.fits')
    bis = append(bis,header_bis['ESO DRS BIS SPAN'])
    fwhm = append(fwhm,header_bis['ESO DRS CCF FWHM'])
    contrast = append(contrast,header_bis['ESO DRS CCF CONTRAST'])
    rv = append(rv,header_bis['ESO DRS CCF RVC'])


keys = ['rv']
format = '%.2f\t'
for i in arange(len(CCF)):
    keys.append('CCF_'+str(i))
    format += '%.1f\t'

dd = dict([(keys[ii],(CCF[ii-1]).tolist()) for ii in arange(1,len(keys),1)])
dd['rv'] = vrad[0]
format = format[:-1]+'\n'
pyrdb.write_rdb('data/'+star+'_CCF.rdb',dd,keys,format)

keys = ['time','CCF','rv','bis_span','fwhm','contrast','sn50']#,'mu_skew','sigma','gamma1']
str_CCF = ['CCF_'+str(i) for i in arange(len(CCF))]
dd = {'time':time,'CCF':str_CCF,'rv':rv,'bis_span':bis,'fwhm':fwhm,'contrast':contrast,'sn50':sn50}#,'mu_skew':mean_skew/1000.,'sigma':sigma_skew,'gamma1':gamma_skew}
format = '%.8f\t%s\t%.5f\t%.5f\t%.5f\t%.5f\t%.2f\n'#\t%.5f\tt%.5f\t%.5f\n'
pyrdb.write_rdb('data/'+star+'_data.rdb',dd,keys,format)

