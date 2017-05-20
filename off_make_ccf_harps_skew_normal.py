#!/usr/bin/env DRSharps
####################################################################
#  off_make_ccf.py [night_directory] [e2ds_fitsfilename] [mask] 
#                  [targetRV] [widthccf] [stepccf] [minorder]
#                  [flux_correction=0/1/template] [wave_fitsfilename]
#
# default  mask=G2
#          targetRV=-99999
#          width=ic_...
#          stepccf=ic_...
#
#  Build New Cross Correlation Function 'ccf'
#   use DRS_DATA_WORKING as output when defined
#
###################################################################
#(@) $Id: off_make_ccf_generic.py,v 3.5 2010/10/13 09:09:13 vltsccm Exp $
###################################################################


import hadmrFITS
import hadgtVISU
import hadmrBIAS
import hadmrDARK
import hadmrEXTOR
import hadmrCDB
import hadmrFLAT
import hadmrLOCOR
import hadrgdCONFIG
import hadmrRV
import fitsio
import fit

import sys,string,time,os,shutil
from numpy.oldnumeric import *
from statis import *
from fitgaus import *
from hadmrMATH import *

import fit_skew_normal


tmp_file="%s"%int(time.time())

execfile(os.getenv('PYTHONSTARTUP'))
execfile(os.getenv('PYTHONSTARTUP2'))

dirfits=os.path.join(dir_data_reduc,arg_night_name)
I_e2ds_fitsfilename=os.path.join(dir_data_reduc,arg_night_name,arg_file_names[0])

if not(locals().has_key('working_dir')):
    working_dir=dirfits


if os.path.exists(working_dir)==0:
   WLOG('error',log_opt,'Redirected outpout directory: '+working_dir+' does not exist')
   exit(1)
else:
   WLOG('info',log_opt,'Using redirected output directory: '+working_dir)


O_e2ds_fitsfilename=os.path.join(working_dir,arg_file_names[0])


WLOG('info',log_opt,'Now processing new CCF on '+I_e2ds_fitsfilename)

if os.path.exists(I_e2ds_fitsfilename):
   insmode=string.strip(hadmrFITS.read_key(I_e2ds_fitsfilename,kw_insmode))
   WLOG('info',log_opt,'INSTRUMENT MODE: '+insmode)
   WLOG('',log_opt,'Reading existing E2DS file '+os.path.split(I_e2ds_fitsfilename)[1])
else:
   WLOG('error',log_opt,'Files  '+I_e2ds_fitsfilename+' missing')
   exit(0)
fiber=''
if (string.find(I_e2ds_fitsfilename,'A.fits') > -1): fiber='A'
if (string.find(I_e2ds_fitsfilename,'B.fits') > -1): fiber='B'
if (fiber==''): 
   WLOG('error',log_opt,'Fiber cannot be identified from filename')
   exit(1)

wave={}
e2dsff={}
blaze={}
param_ll={}

#
# Read image keywords
#
filename=I_e2ds_fitsfilename
execfile(ic_exec_kw_list)

try:
   _kw_out_=fitsio.read_keys(filename,kw)
except IOError:
   WLOG('error',log_opt,'I/O Error on file: '+filename)
   exit(1)

WLOG('info',log_opt,'Now processing Image TYPE: '+_kw_out_[kw_dprtype][0]+'  with '+process_running+' recipe')
dprtype=_kw_out_[kw_dprtype][0]

execfile(ic_exec_kw_allocate)
execfile(ic_exec_kw_display)


e2dsff[fiber],nx,ny=hadmrFITS.read_data(I_e2ds_fitsfilename)


blaze_file=os.path.join(dirfits,hadmrFITS.read_key(I_e2ds_fitsfilename,kw_BLAZE_FILE[0]))
if os.path.exists(blaze_file): 
    blaze[fiber],fx,fy=hadmrFITS.read_data(blaze_file)
else:
    WLOG('error',log_opt,'blaze file: '+blaze_file+' does not exist, FATAL!')
    exit(0)

kw=[kw_tarrv,kw_BERV[0],kw_BJD[0],kw_BERVMX[0],kw_snref,kw_DRIFT_VR_USED[0],kw_CCD_SIGDET[0]]
_kw_out_=fitsio.read_keys(I_e2ds_fitsfilename,kw,hdu=-1,verbose=0)
berv=_kw_out_[kw_BERV[0]][0]
berv_max=_kw_out_[kw_BERVMX[0]][0]
bjd=_kw_out_[kw_BJD[0]][0]
WLOG('info',log_opt,'Barycentric Earth RV correction: %9.5f km/s' %(berv))
kw_BERV[1]=berv
kw_BJD[1]=bjd
kw_BERVMX[1]=berv_max    
snref=_kw_out_[kw_snref][0]
drift_instr=float(_kw_out_[kw_DRIFT_VR_USED[0]][0])
ccdsigdet=_kw_out_[kw_CCD_SIGDET[0]][0]
noise=ccdsigdet*sqrt(ic_extmeanzone)


try:
   wave_fitsfilename = arg_file_names[7]
except IndexError:
   wave[fiber],param_ll[fiber],paramx=hadmrFITS.get_e2ds_ll(I_e2ds_fitsfilename,_kw_TH_CAL_)

if locals().has_key('wave_fitsfilename'):
   WLOG('info',log_opt,'Reading wavelength calibration from file: '+wave_fitsfilename)
   dic_db = {}
   dic_db['TH_A'] = wave_fitsfilename
   dic_db['THREF_A'] = string.replace(wave_fitsfilename,'wave_A','e2ds_A')
   dic_db['THCCF_A'] = string.replace(wave_fitsfilename,'wave_A','ccf_TH_A')
   dic_db['THRV_A'] = string.replace(wave_fitsfilename,'wave_A.fits','ccf_TH_A.tbl')
   dic_db['TH_B'] = string.replace(wave_fitsfilename,'wave_A','wave_B')
   dic_db['THREF_B'] = string.replace(wave_fitsfilename,'wave_A','e2ds_B')
   dic_db['THCCF_B'] = string.replace(wave_fitsfilename,'wave_A','ccf_TH_B')
   dic_db['THRV_B'] = string.replace(wave_fitsfilename,'wave_A.fits','ccf_TH_B.tbl')
   test_global = 1
   test_local = 1
   for key in dic_db.keys():
      test_global *= os.path.exists(os.path.join(dir_data_reduc,arg_night_name,dic_db[key]))
      test_local *= os.path.exists(os.path.join(working_dir,dic_db[key]))
   if test_global == 1:
      wavedir = os.path.join(dir_data_reduc,arg_night_name)
   elif test_local == 1:
      wavedir = working_dir
   else:
      WLOG('error',log_opt,'Impossible to find wavelength calibration files')
      exit(1)
   WLOG('',log_opt,'Doing instrumental drift computation:') 
   B_e2ds_fitsfilename = string.replace(I_e2ds_fitsfilename,'e2ds_A','e2ds_B')
   try:
      spe,nx,ny=hadmrFITS.read_data(B_e2ds_fitsfilename)
   except IOError:
      WLOG('error',log_opt,'Impossible to find simultaneous reference spectrum')
      exit(1)
   th_file=os.path.join(wavedir,dic_db['THREF_A'])
   thref_file=os.path.join(wavedir,dic_db['THREF_B'])
   thccfref_file=os.path.join(wavedir,dic_db['THCCF_B'])
   thRVref_file=os.path.join(wavedir,dic_db['THRV_B'])
   wave['B'],param_ll['B'],paramx=hadmrFITS.get_e2ds_ll(thref_file,_kw_TH_CAL_)
   wave[fiber],param_ll[fiber],paramx=hadmrFITS.get_e2ds_ll(th_file,_kw_TH_CAL_)
   th_RVref=hadmrRV.get_RVCCF(thRVref_file)
   WLOG('',log_opt,'Reading E2DS reference spectrum '+thref_file)
   speref,nx,ny=hadmrFITS.read_data(thref_file)
   noise=ccdsigdet*sqrt(ic_extmeanzone)
   
   if ic_drift_algo == 'Bouchy':
   
      WLOG('info',log_opt,'Bouchy algorithm is used for drift estimate')
      dvrmsref2,meanpondref2=hadmrRV.deltavrms2d(speref,wave['B'],noise)
      dvrmsspe2,meanpondspe2=hadmrRV.deltavrms2d(spe,wave['B'],noise)
      WLOG('info',log_opt,'Photon noise uncertainty on Thorium = %.3f m/s' %(meanpondspe2))
      
      WLOG('',log_opt,'Normalizing spectrum and doing cosmic correction ')
      spen,u,cpt=hadmrRV.renorm_cosmics2d(speref,spe)
      
      WLOG('',log_opt,'Computing drift')
      vr=hadmrRV.calcvrdrifts2d(speref,spen,wave['B'],noise)
      
      rap=sum(u)/len(u)
      wref=1.0/dvrmsref2
      
      meanvr=sum(vr*wref)/sum(wref)
      err_meanvr=sqrt(dvrmsref2+dvrmsspe2+ic_disp_drift**2)
      #
      # sig-clip
      #
      sig_val=abs(vr-meanvr)/err_meanvr
      #hadgtVISU.megaplot(arange(len(vr)),vr,ymin=-50,ymax=50)
      while sum(greater(sig_val,4)):
   	   wref[argsort(sig_val)[-1]]=0
   	   meanvr=sum(vr*wref*greater(wref,0))/sum(wref)
   	   sig_val=abs(vr-meanvr)*greater(wref,0)/err_meanvr
      nborderout=len(vr)-sum(greater(wref,0))
      
      WLOG('info',log_opt,'On thorium  Drift=%.2f m/s,  %i cosmics, Flux ratio=%.3f, %s rejected orders ' %(-1*meanvr,cpt,rap,nborderout))
      wref=1.0/dvrmsref2
      
      #
      # QC on drift
      #
      if rap>qc_drift_rapflu or rap<(1./qc_drift_rapflu): 
         BAD_DRIFT='Flux ratio between ref. and current thorium is inadequate for bouchy drift' 
      if abs(meanvr)>qc_drift_max:
         BAD_DRIFT='Spectro drift [%5.1f] (m/s) is too large for bouchy technique'%(meanvr)
      if nborderout>qc_drift_nborderout: 
         BAD_DRIFT='Too many orders discarded for bouchy drift measurement'
      if meanpondspe2>qc_max_err_on_simulthorium:
         BAD_DRIFT='RV Photon noise error on Thorium flux too high for drift computation'
      
      #
      # Fill in drift keywords
      #
      kw_DRIFT_REF[1]=dic_db['THREF_B']
      kw_DRIFT_REF_CCF[1]=dic_db['THCCF_B']
      kw_DRIFT_REF_RV[1]=0.
      kw_DRIFT_VR[1]=int(1000.*(-1)*meanvr)/1000.
      kw_DRIFT_VR_CCF[1]=0.
      kw_DRIFT_NBCOS[1]=cpt
      kw_DRIFT_RFLUX[1]=int(1000.*rap)/1000.
      kw_DRIFT_NBORDKILL[1]=nborderout
      
      if locals().has_key('BAD_DRIFT'):
         WLOG('error',log_opt,'Drift correction cannot be applied: '+BAD_DRIFT)
         kw_DRIFT_VR_USED[1]=0.
         kw_DRIFT_ALGO[1]='NONE'
         kw_DRIFT_NOISE[1]=ic_typical_drift
         kw_DRIFT_QC[1]='FAILED'
      else:
         kw_DRIFT_VR_USED[1]=kw_DRIFT_VR[1]
         kw_DRIFT_ALGO[1]='Bouchy'
         kw_DRIFT_NOISE[1]=int(1000.*meanpondspe2)/1000.
         kw_DRIFT_QC[1]='PASSED'
      
   else:
      
      WLOG('info',log_opt,'CCF algorithm is used for drift estimate')
      dvrmsref2,meanpondref=hadmrRV.deltavrms2d(speref,wave['B'],noise)
      wref=1.0/dvrmsref2
      
      cor_fitsfilename={\
          'B': os.path.join(working_dir,string.replace(arg_file_names[0],\
                        '_e2ds_A','_ccf_TH_B'))}   
      #th simult
      RV_template=os.path.join(dir_drs_config,ic_thorium_mask+'.mas')
      ll_mask_D,ll_mask_ctr,w_mask=hadmrRV.get_mask(RV_template,0,1,log_opt)
      RV_CCF=arange(-ic_thorium_CCF_half_width,ic_thorium_CCF_half_width+ic_thorium_CCF_step,ic_thorium_CCF_step)
      CCF,CCF_max,pix_passed_all,tot_line,ll_range_all,CCF_noise=\
        hadmrRV.coravelation(wave['B'],spe,param_ll['B'],spe*0.+1.,RV_CCF,ll_mask_D,ll_mask_ctr,w_mask,0,0,noise,'.th',1,ic_debug)
      #ref th
      RV_CCF_ref=arange(-ic_thorium_CCF_half_width,ic_thorium_CCF_half_width+ic_thorium_CCF_step,ic_thorium_CCF_step)
      CCF_ref,CCF_max_ref,pix_passed_all_ref,tot_line_ref,ll_range_all_ref,CCF_noise=\
        hadmrRV.coravelation(wave['B'],speref,param_ll['B'],speref*0.+1.,\
                          RV_CCF_ref,ll_mask_D,ll_mask_ctr,w_mask,0,0,noise,'.thref',1,ic_debug)
      #
      # reject bad order by sig-clip
      #
      RV_o=[]
      RV_ref=[]
      for o in (range(len(CCF))):
        CCF_res_o,CCF_fit_o=hadmrRV.fit_CCF(RV_CCF,CCF[o],1)
        CCF_res_ref,CCF_fit_ref=hadmrRV.fit_CCF(RV_CCF_ref,CCF_ref[o],1)
        RV_o.append(CCF_res_o[1])
        RV_ref.append(CCF_res_ref[1])
      flag=1.0+0*wref
      Diff_rv=(RV_o*flag-RV_ref*flag)
      #wref=th_RVref[:,1]
      mean=sum(Diff_rv*wref)/sum(wref)
      mean2=sum(Diff_rv**2*wref)/sum(wref)
      sig=sqrt(mean2-mean**2)
      sig_val=abs((Diff_rv-mean)/sig)
      #hadgtVISU.megaplot(arange(len(Diff_rv)),Diff_rv,ymin=mean-4*sig,ymax=mean+4*sig)
      #hadgtVISU.megaplot(arange(len(Diff_rv)),RV_o,ymin=-0.050,ymax=0.050)
      while sum(greater(sig_val,4.)*flag):
        flag[argsort(sig_val*flag)[-1]]=0
        mean=sum(Diff_rv*wref*flag)/sum(wref*flag)
        mean2=sum(Diff_rv**2*wref*flag)/sum(wref*flag)
        sig=sqrt(mean2-mean**2)
        sig_val=abs((Diff_rv-mean)/sig)
      for o in (range(len(CCF))):
         CCF[o]=CCF[o]*flag[o]
      nborder_rm=len(CCF)-sum(flag)
      average_CCF=sum(CCF)
      normalized_CCF=average_CCF/max(average_CCF)
      CCF_res,CCF_fit=hadmrRV.fit_CCF(RV_CCF,normalized_CCF,1)
      # redo centered   
      RV_CCF=arange(CCF_res[1]-ic_thorium_CCF_half_width,CCF_res[1]+ic_thorium_CCF_half_width+ic_thorium_CCF_step,ic_thorium_CCF_step)
      CCF,CCF_max,pix_passed_all,tot_line,ll_range_all,CCF_noise=\
        hadmrRV.coravelation(wave['B'],spe,param_ll['B'],spe*0.+1.,\
                          RV_CCF,ll_mask_D,ll_mask_ctr,w_mask,0,0,noise,string.replace(cor_fitsfilename['B'],'.fits','.tbl'),1,ic_debug)
      for o in (range(len(CCF))):
         CCF[o]=CCF[o]*flag[o]
      
      average_CCF=sum(CCF)
      normalized_CCF=average_CCF/max(average_CCF)
      CCF_res,CCF_fit=hadmrRV.fit_CCF(RV_CCF,normalized_CCF,1)
      CCF_res[0]=CCF_res[0]/(1.+CCF_res[3])
      WLOG('',log_opt,'Correlation on simult thorium: '+\
   	       ' C='+"%.1f"%(abs(100*CCF_res[0]))+'[%]'+\
   	       ' RV='+"%.5f"%(CCF_res[1])+'[km/s]'+\
   	       ' FWHM='+"%.4f"%(CCF_res[2]*2.3548)+'[km/s]'+\
   	       ' maxcpp='+"%.1f"%(sum(CCF_max)/sum(pix_passed_all)))
      
      err=sqrt(average_CCF+sum(pix_passed_all)*noise**2)
      result=fit.gauss(RV_CCF,average_CCF,err)
      
      #
      # rebuild ref RV
      #
      RV_CCF_ref=arange(-ic_thorium_CCF_half_width,ic_thorium_CCF_half_width+ic_thorium_CCF_step,ic_thorium_CCF_step)
      CCF_ref,CCF_max_ref,pix_passed_all_ref,tot_line_ref,ll_range_all_ref,CCF_noise=\
        hadmrRV.coravelation(wave['B'],speref,param_ll['B'],speref*0.+1.,\
                          RV_CCF_ref,ll_mask_D,ll_mask_ctr,w_mask,0,0,noise,'.tmp',1,ic_debug)
      for o in (range(len(CCF))):
         CCF_ref[o]=CCF_ref[o]*flag[o]
      
      average_CCF=sum(CCF_ref)
      normalized_CCF_ref=average_CCF/max(average_CCF)
      CCF_res_ref,CCF_fit_ref=hadmrRV.fit_CCF(RV_CCF_ref,normalized_CCF_ref,1)
      CCF_res_ref[0]=CCF_res_ref[0]/(1.+CCF_res_ref[3])
      WLOG('',log_opt,'Correlation on ref thorium : '+\
   	       ' C='+"%.1f"%(abs(100*CCF_res_ref[0]))+'[%]'+\
   	       ' RV='+"%.5f"%(CCF_res_ref[1])+'[km/s]'+\
   	       ' FWHM='+"%.4f"%(CCF_res_ref[2]*2.3548)+'[km/s]'+\
   	       ' maxcpp='+"%.1f"%(sum(CCF_max_ref)/sum(pix_passed_all_ref)))
      
      err=sqrt(average_CCF+sum(pix_passed_all_ref)*noise**2)
      result_ref=fit.gauss(RV_CCF_ref,average_CCF,err)
      
      meanvr2=(CCF_res[1]-CCF_res_ref[1])*1000
      WLOG('info',log_opt,'On thorium CCF Drift =%.2f m/s  , %d  order(s) removed'  %(meanvr2,nborder_rm))
      summary='Drift on CCF= %.3f m/s'%(meanvr2)
      
      #
      # archive in fits
      #                
      WLOG('',log_opt,'Archiving CCF on file: '+cor_fitsfilename['B'])
      hadmrFITS.Save_ccf(RV_CCF,CCF,sum(CCF),cor_fitsfilename['B'],ic_thorium_mask,\
   	   sum(CCF_max)/sum(pix_passed_all),sum(tot_line),abs(100*CCF_res[0]),CCF_res[1],CCF_res[2]*2.3548,CCF_res[1],_kw_CCF)
      
      #
      # QC on drift
      #
      if abs(100*CCF_res[0])<qc_drift_CCF_minC:
         BAD_DRIFT='Low contrast of thorium CCF, check thorium flux, drift value is uncertain'
      if abs(meanvr2)>qc_drift_max:
         WLOG('warning',log_opt,'Drift is surprisingly large ('+"%.0f"%meanvr2+' [m/s]), I would recommend to do a TH calib as soon as possible')
      
      #
      # Fill in drift keywords
      #
      kw_DRIFT_REF[1]=dic_db['THREF_B'][1]
      kw_DRIFT_REF_CCF[1]=dic_db['THCCF_B'][1]
      kw_DRIFT_REF_RV[1]=int(1000.*CCF_res_ref[1])/1000.
      kw_DRIFT_VR[1]=0.
      kw_DRIFT_VR_CCF[1]=int(1000.*meanvr2)/1000.
      kw_DRIFT_NBCOS[1]=0
      kw_DRIFT_RFLUX[1]=0.
      kw_DRIFT_NBORDKILL[1]=nborder_rm
      
      if locals().has_key('BAD_DRIFT'):
         WLOG('error',log_opt,'Drift can not be computed: '+BAD_DRIFT)
         kw_DRIFT_VR_USED[1]=0.
         kw_DRIFT_ALGO[1]='NONE'
         kw_DRIFT_NOISE[1]=ic_typical_drift
         kw_DRIFT_QC[1]='FAILED'
      else:
         kw_DRIFT_VR_USED[1]=kw_DRIFT_VR_CCF[1]
         kw_DRIFT_ALGO[1]='CCF'
         kw_DRIFT_NOISE[1]=sqrt(result[7]**2+result_ref[7]**2)*1000.
         kw_DRIFT_QC[1]='PASSED'
      
   drift_instr = kw_DRIFT_VR_USED[1]



try:
   ccf_mask=arg_file_names[1]
   if ccf_mask == 'kw':
      ccf_mask=hadmrRV.get_mask_name(target_Sptype,ic_default_mask,os.path.join(dir_drs_config,ic_CCF_template_db),log_opt)
except IndexError:
   ccf_mask=hadmrRV.get_mask_name(target_Sptype,ic_default_mask,os.path.join(dir_drs_config,ic_CCF_template_db),log_opt)

try:
   target_RV=arg_file_names[2]
   if target_RV == 'kw': target_RV=float(_kw_out_[kw_tarrv][0])
   else: target_RV=float(arg_file_names[2])
except IndexError:
   target_RV=float(_kw_out_[kw_tarrv][0])

############
# xavier dumusque correction
# to prevent oversampling the CCF
#
#try:
#   ccf_width=float(arg_file_names[3])
#except IndexError:
#   ccf_width=ic_ccf_width
ccf_width = ic_plate_scale * 25.
############


############
# xavier dumusque correction
# to prevent oversampling the CCF
#try:
#   ccf_step=float(arg_file_names[4])
#except IndexError:
#   ccf_step=ic_ccf_step
ccf_step = ic_plate_scale
############

try:
   ccf_order_range=arange(int(arg_file_names[5]),ic_ccf_order_range[fiber][-1]+1)
except IndexError:
   ccf_order_range=ic_ccf_order_range[fiber]

try:
   do_flux_correction=int(arg_file_names[6])
except IndexError:
   do_flux_correction=ic_do_flux_correction
except ValueError:
   do_flux_correction=1
   target_Sptype=arg_file_names[6]

#
# Cosmic filter correction
#

newcosmic=0
WLOG('',log_opt,'Doing Cosmic filtering')
e2dsff[fiber]=e2dsff[fiber]*1.0
for no in range(0,len(e2dsff[fiber])):
    e2dsff[fiber][no],nbcos2=hadmrEXTOR.cosmic_filter(e2dsff[fiber][no]*1.0)
    newcosmic=newcosmic+nbcos2
    
WLOG('info',log_opt,' %i cosmic hit corrected'%(newcosmic))


##################
# DO CORRELATION #
##################

#
# get RV template
#
WLOG('info',log_opt,'Template used for CCF computation: '+ccf_mask)
RV_template=os.path.join(dir_drs_config,ccf_mask+'.mas')
ll_mask_D,ll_mask_ctr,w_mask=hadmrRV.get_mask(RV_template,ic_w_mask_min,ic_mask_width,log_opt)

#
# Relative flux correction
#

if do_flux_correction:
   flux_template = hadmrRV.get_flux_template(target_Sptype,insmode,fiber,os.path.join(dir_drs_config,ic_flux_template_db),log_opt)
   try:
      flux_ref = hadmrFITS.read_rdb(os.path.join(dir_drs_config,flux_template))['flux']
   except IOError:
      WLOG('info',log_opt,'Flux correction not performed: no flux template for spectral type '+target_Sptype)
      do_flux_correction = 0

corr = wave[fiber]*0.+1.
p_corr = array([1.,0.,0.,0.,0.,0.])

if do_flux_correction:
   
   nx = len(e2dsff[fiber][0])
   ny = len(e2dsff[fiber])
   flux = zeros(ny,'d')
   flux_ref = asarray(flux_ref)
   central_ll = zeros(ny,'d')
   
   for i in range(ny):
      flux[i] = sum(sort(e2dsff[fiber][i,nx/4:3*nx/4])[:-int(nx/400)])
      central_ll[i] = wave[fiber][i,nx/2]
   
   sn_min = 100.
   qq = greater(flux,0.5*sn_min**2*(1.+sqrt(1.+4*nx/2*noise**2/sn_min**2)))
   
   if sum(qq) > (0.8*ny):
      flux = compress(qq,flux)
      flux_ref = compress(qq,flux_ref)
      central_ll = compress(qq,central_ll)
      corr_factor = flux/flux_ref/stats.mean(flux/flux_ref)
      p_corr,sig_p,mod,chisq = fit.poly(central_ll,corr_factor,ones(len(flux),'d'),6)
      corr = wave[fiber]*0.
      for i in range(len(p_corr)): corr += p_corr[i]*wave[fiber]**i
      if (corr.min() > 0.25) & (corr.max() < 3):
         e2dsff[fiber] = e2dsff[fiber]/corr
         WLOG('info',log_opt,'Relative flux correction performed with min/max weights = '+\
             "%1.2f"%corr.min()+'/'+"%1.2f"%corr.max())
      else:
         WLOG('warning',log_opt,'Flux correction not performed: mismatch with flux template too large')
   else:
      WLOG('info',log_opt,'Global flux too low to perform relative flux correction')


cor_fitsfilename={\
     'A': os.path.join(working_dir,string.replace(arg_file_names[0],\
                       '_e2ds','_ccf_'+ccf_mask)),\
     'B': os.path.join(working_dir,string.replace(arg_file_names[0],\
                       '_e2ds','_ccf_'+ccf_mask))}

res_tbl_file=string.replace(cor_fitsfilename[fiber],'.fits','.tbl')


if int(target_RV)==-99999:
   WLOG('',log_opt,'Searching CCF is ongoing')
   RV_CCF=arange(-ic_ccf_Wide_width,ic_ccf_Wide_width+ic_ccf_Wide_step,ic_ccf_Wide_step)
   CCF,CCF_max,pix_passed_all,tot_line,ll_range_all,CCF_noise=\
       hadmrRV.coravelation(wave[fiber],e2dsff[fiber],param_ll[fiber],blaze[fiber],RV_CCF,ll_mask_D,\
            ll_mask_ctr,w_mask,berv,berv_max,noise,res_tbl_file,0,ic_debug)
   target_RV=RV_CCF[argsort(sum(CCF))[0]]
   WLOG('',log_opt,'Guessed RV is at '+"%6.1f"%target_RV+'[km/s]')

WLOG('',log_opt,'Computing CCF at RV= '+"%6.1f"%target_RV+'[km/s]')

RV_CCF=arange(target_RV-ccf_width,target_RV+ccf_width+ccf_step,ccf_step)

if (ic_do_blaze_correct):
   WLOG('',log_opt,'Doing Blaze correction')
   effblaze=blaze[fiber]
else:
   WLOG('',log_opt,'Blaze correction is disabled')
   effblaze=blaze[fiber]*0.+1.

CCF,CCF_max,pix_passed_all,tot_line,ll_range_all,CCF_noise=\
       hadmrRV.coravelation(wave[fiber],e2dsff[fiber],param_ll[fiber],effblaze,RV_CCF,ll_mask_D,\
            ll_mask_ctr,w_mask,berv,berv_max,noise,res_tbl_file,0,ic_debug)

average_CCF=sum(CCF[ccf_order_range])
normalized_CCF=average_CCF/max(average_CCF)
CCF_res,CCF_fit=hadmrRV.fit_CCF(RV_CCF,normalized_CCF,0)
CCF_res[0]=CCF_res[0]/(1.+CCF_res[3])
maxcpp=sum(CCF_max)/sum(pix_passed_all)

WLOG('info',log_opt,'Correlation Normal: '+\
	    ' C='+"%.1f"%(abs(100*CCF_res[0]))+'[%]'+\
	    ' RV='+"%.5f"%(CCF_res[1])+'[km/s]'+\
	    ' FWHM='+"%.4f"%(CCF_res[2]*2.3548)+'[km/s]'+\
	    ' maxcpp='+"%.1f"%(sum(CCF_max)/sum(pix_passed_all)))

############
# xavier dumusque correction
# to fit the skew normal distribution to the data
mod_skew,solution,init_mod_skew,chi2,chi2_red,e_skew,w_skew,a_skew = fit_skew_normal.fit_skew_distribution_for_CCF(RV_CCF,normalized_CCF,CCF_noise,3000.,CCF_res[1],CCF_res[2]*2.3548,method='leastsq')

WLOG('info',log_opt,'Correlation Skew Normal: '+\
	    ' C='+"%.1f"%(abs(100*CCF_res[0]))+'[%]'+\
	    ' RV='+"%.5f"%(CCF_res[1])+'[km/s]'+\
	    ' FWHM='+"%.4f"%(CCF_res[2]*2.3548)+'[km/s]'+\
	    ' maxcpp='+"%.1f"%(sum(CCF_max)/sum(pix_passed_all)))
RV_skew=CCF_res[1]-drift_instr/1000.     
############

RV=CCF_res[1]-drift_instr/1000.
WLOG('info',log_opt,'RV corrected from drift (%.2f m/s) is: %.5f [km/s]'%(drift_instr,RV))


#
# Photon noise uncertainty on fiber A
#
dvrmsref2,meanpondref=hadmrRV.deltavrms2d(e2dsff[fiber],wave[fiber],noise)
if locals().has_key('ic_exec_correct_meanpondref'): 
    WLOG('',log_opt,'Correcting error estimate')
    execfile(ic_exec_correct_meanpondref)

CCF_noise_tot=sqrt(sum(CCF_noise[ccf_order_range]**2))
CCF_slope=(average_CCF[2:]-average_CCF[:-2])/(RV_CCF[2:]-RV_CCF[:-2])
ccf_oversampling=ic_plate_scale/ccf_step
indexlist=map(int,arange(round(len(CCF_slope)/ccf_oversampling))*ccf_oversampling)
qq=zeros(len(CCF_slope))
for i in indexlist: qq[i]=1
ccf_noise=(sum(compress(qq,CCF_slope**2)/compress(qq,CCF_noise_tot[1:-1]**2)))**(-0.5)

WLOG('info',log_opt,'Estimated RV accuracy on stellar spectrum: %.2f[m/s]' %(meanpondref))
WLOG('info',log_opt,'Estimated RV accuracy on CCF             : %.2f[m/s]' %(ccf_noise*1000.))

#
# archive in fits
#
WLOG('',log_opt,'Archiving CCF on file: '+cor_fitsfilename[fiber])
hadmrFITS.Save_ccf(RV_CCF,CCF,average_CCF,cor_fitsfilename[fiber],ccf_mask,\
sum(CCF_max)/sum(pix_passed_all),sum(tot_line[ccf_order_range]),abs(100*CCF_res[0]),CCF_res[1],CCF_res[2]*2.3548,RV,_kw_CCF)

# Add CCF photon noise keyword
kw_CCF_NOISE[1]=ccf_noise
hadmrFITS.write_newkey(cor_fitsfilename[fiber],kw_CCF_NOISE)

#
# copy e2ds file kw into CCF
#
hadmrFITS.copy_key(I_e2ds_fitsfilename,cor_fitsfilename[fiber])

# Update DVRMS keyword
kw_DVRMS[1]=meanpondref
hadmrFITS.update_key(cor_fitsfilename[fiber],kw_DVRMS)

# Add flux correction keywords
kw_FLUX_CORR_MIN[1]=corr.min()
kw_FLUX_CORR_MAX[1]=corr.max()
hadmrFITS.write_newkey(cor_fitsfilename[fiber],kw_FLUX_CORR_MIN)
hadmrFITS.write_newkey(cor_fitsfilename[fiber],kw_FLUX_CORR_MAX)
for i in range(len(p_corr)):
   hadmrFITS.write_newkey(cor_fitsfilename[fiber],[kw_FLUX_CORR_COEFF[0]+str(i),p_corr[i],kw_FLUX_CORR_COEFF[2]])

#
# Compute bisector
#
depth,bis,span,v0 = hadmrRV.compute_bis(RV_CCF,average_CCF)
bis_fitsfilename=string.replace(cor_fitsfilename[fiber],'ccf','bis')
WLOG('',log_opt,'Archiving bisector on file: '+bis_fitsfilename)
hadmrFITS.Save_bis(depth,bis,span,v0,bis_fitsfilename,_kw_BIS)
hadmrFITS.copy_key(cor_fitsfilename[fiber],bis_fitsfilename)

# Update drift keywords if necessary
if locals().has_key('wave_fitsfilename'):
   hadmrFITS.update_key(cor_fitsfilename[fiber],kw_DRIFT_REF)
   hadmrFITS.update_key(cor_fitsfilename[fiber],kw_DRIFT_VR)
   hadmrFITS.update_key(cor_fitsfilename[fiber],kw_DRIFT_NBCOS)
   hadmrFITS.update_key(cor_fitsfilename[fiber],kw_DRIFT_RFLUX)
   hadmrFITS.update_key(cor_fitsfilename[fiber],kw_DRIFT_NBORDKILL)
   hadmrFITS.update_key(cor_fitsfilename[fiber],kw_DRIFT_NOISE)
   hadmrFITS.update_key(cor_fitsfilename[fiber],kw_DRIFT_REF_CCF)
   hadmrFITS.update_key(cor_fitsfilename[fiber],kw_DRIFT_REF_RV)
   hadmrFITS.update_key(cor_fitsfilename[fiber],kw_DRIFT_VR_CCF)
   hadmrFITS.update_key(cor_fitsfilename[fiber],kw_DRIFT_VR_USED)
   hadmrFITS.update_key(cor_fitsfilename[fiber],kw_DRIFT_ALGO)
   hadmrFITS.update_key(cor_fitsfilename[fiber],kw_DRIFT_QC)


############
# display
############

if os.environ.has_key('DRS_PLOT'):  
   summary='C='+"%.1f"%(abs(100*CCF_res[0]))+'%'+\
	    '  RVC='+"%.5f"%(RV)+' km/s'+\
	    '  FWHM='+"%.4f"%(CCF_res[2]*2.3548)+' km/s'

if os.environ.has_key('DRS_PLOT'):  
   if os.environ['DRS_PLOT']=='trigger':
      hadgtVISU.printfic('/tmp/plotA.'+tmp_file,RV_CCF,normalized_CCF)
      WLOG('graph',log_opt,"hadgtVISU.megaplotfic('/tmp/plotA."+tmp_file+"',xmin="+'%s'%int(RV_CCF[0]-1)+",xmax="+'%s'%int(RV_CCF[-1]+1)+",ymin=0.25,ymax=1.05)")
      hadgtVISU.printfic('/tmp/plotB.'+tmp_file,RV_CCF,CCF_fit)
      WLOG('graph',log_opt,"hadgtVISU.megaplotfic('/tmp/plotB."+tmp_file+"',xmin="+'%s'%int(RV_CCF[0]-1)+",xmax="+'%s'%int(RV_CCF[-1]+1)+",ymin=0.25,ymax=1.05,overplot=1,xtitle='km/s',ytitle='CCF',title='Average CCF normalized of "+target_name+" with TPL= "+ccf_mask+"\\\\n "+summary+"')")     
   else:
      hadgtVISU.megaplot(RV_CCF,normalized_CCF,xmin=RV_CCF[0]-1,xmax=RV_CCF[-1]+1,ymin=0.25,ymax=1.05)
      hadgtVISU.megaplot(RV_CCF,CCF_fit,xmin=RV_CCF[0]-1,xmax=RV_CCF[-1]+1,ymin=0.25,ymax=1.05,overplot=1,xtitle='km/s',ytitle='CCF',title='Average CCF normalized of '+target_name+' with TPL= '+ccf_mask+'\\n '+summary)


#
# do tar archive for NGACS system
#
if  os.environ.has_key('HAR_NGAS_DIR') :
   tar_in_files=os.path.join(dir_data_reduc,arg_night_name,string.split(arg_file_names[0],'_')[0]+'*')
   failed=os.system('tar cvf ./reduced_tar.tmp '+tar_in_files+' >& /tmp/.drs')
   if failed:
      WLOG('warning',log_opt,'Unable to write tar file ./reduced_tar.tmp for NGAS archive')
   else:
      WLOG('',log_opt,'tmp tar file ./reduced_tar.tmp done')
   
   if (os.path.exists('./reduced_tar.tmp')):
      try:
         shutil.copy('./reduced_tar.tmp',os.path.join(os.environ['HAR_NGAS_DIR'],string.split(arg_file_names[0],'_')[0]+'.tar'))
      except OSError:
         WLOG('warning',log_opt,'I/O error for NGAS archive')
      else:
         WLOG('',log_opt,'tar file moved to '+os.environ['HAR_NGAS_DIR'])
else:
   WLOG('',log_opt,'NGAS archiving system is disable')


WLOG('info',log_opt,'Recipe '+process_running+' is terminated')


