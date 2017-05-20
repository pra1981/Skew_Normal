#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib.pyplot as pl


from pylab import *
from scipy.special import erf
from scipy.stats import norm
from scipy.optimize import leastsq

from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm

import sys
sys.path.append('/Users/xav/uni/phd/program/functions/')

import stats
import functions
import fit2
import pymc
import random



# For more information about the skew distribution, see "http://azzalini.stat.unipd.it/SN/Intro/intro.html"
#and have a look at "test/skew_distribution_for_CCF"
#def pdf(x):
#    return 1/sqrt(2*pi) * exp(-x**2/2.)
#
#def cdf(x):
#    return (1 + erf(x/sqrt(2))) / 2.
#
#def skew(x,e=0,w=1,a=0):
#    t = (x-e) / (1.*w)
#    return 2./w * pdf(t) * cdf(a*t)
#    #return 2./w * norm.pdf(t) * norm.cdf(a*t)
### You can of course use the scipy.stats.norm versions
### return 2 * norm.pdf(t) * norm.cdf(a*t)


n = 2**10

e = 5.0 # location
w = 5000. # scale

x = linspace(-50000,50000+1e-10,n)

#skew_vect = arange(-1.,1.+1e-10,0.01)
skew_vect = arange(-0.5,0.5+1e-10,0.001)
mean_p,mean_p_corr,mean_fit_p,mu_fitted,gamma_fitted,vrad_norm_vect,span_vect,a_fitted = [],[],[],[],[],[],[],[]
mu_real,sigma_real,gamma_real = [],[],[]
for i,a in enumerate(skew_vect):

    mu,sigma,gamma,gamma_2 = functions.transform_skew_para_from_dp_to_cp(e,w,a)
    mu_real = append(mu_real,mu)
    sigma_real = append(sigma_real,sigma)
    gamma_real = append(gamma_real,gamma)
    
    correct_mean_offset = functions.correct_mean_offset(w,a)
    
    if i%1000 == 0:
        print i,' over ',len(skew_vect)
    cte = 0.
    ampli = -3000.
    p = cte-ampli*functions.skew(x,e,w,a)
    error = array([random.gauss(0,0.0002) for j in arange(len(x))])
    p += error
    p_cp = cte-ampli*functions.skew_cp(x,mu,sigma,gamma) + error
    
    mean_p = append(mean_p,stats.wmean(x,cte-p))
    p_corr = cte-ampli*functions.skew(x,e-correct_mean_offset,w,a)
    mean_p_corr = append(mean_p_corr,stats.wmean(x,cte-p_corr))

    mod,solution,init_mod,chi2,chi2_red,e_fit,w_fit,a_fit = functions.fit_skew_distribution_for_CCF(x,p,-ampli,e,2*sqrt(2.*log(2.))*w)
#    mean_fit_p = append(mean_fit_p,stats.wmean(x,solution[0]-mod))
#    mu_fitted = append(mu_fitted,solution[2])
#    gamma_fitted = append(gamma_fitted,solution[4])
#    mean_fit_p = append(mean_fit_p,stats.wmean(x,solution['cte'].value-mod))
#    mu_fitted = append(mu_fitted,solution['mu'].value)
#    gamma_fitted = append(gamma_fitted,solution['gamma'].value)
    mean_fit_p = append(mean_fit_p,stats.wmean(x,-mod))
    mu_fitted = append(mu_fitted,solution[1])
    gamma_fitted = append(gamma_fitted,solution[3])
    a_fitted = append(a_fitted,a_fit)

#    mod_norm,continuum,contrast,span,vrad_norm,fwhm,depth,bis = functions.compute_bis(x,p/cte)
#    
#    vrad_norm_vect = append(vrad_norm_vect,vrad_norm)
#    span_vect = append(span_vect,span)

    subplot(211)
    plot(x,p)
    plot(x,init_mod,'k')
    plot(x,mod,c='r',ls='--')
    #    plot(x,cte*mod_norm,c='g',ls='-')
    subplot(212)
    plot(x,p_corr)


#skew_vect = arange(-0.14,0.14+1e-10,0.001)
##skew_vect = arange(0.5,0.6,1.)
#mean_p,mean_p_corr,mean_fit_p,mu_fitted,gamma_fitted = [],[],[],[],[]
#mu_real,sigma_real,gamma_real = [],[],[]
#for i,a in enumerate(skew_vect):
#    
#    mu,sigma,gamma,gamma_2 = functions.transform_skew_para_from_dp_to_cp(e,w,a)
#    mu_real = append(mu_real,mu)
#    sigma_real = append(sigma_real,sigma)
#    gamma_real = append(gamma_real,gamma)
#    
#    correct_mean_offset = functions.correct_mean_offset(w,a)
#    
#    if i%1000 == 0:
#        print i,' over ',len(skew_vect)
#    p = functions.skew(x,e,w,a)
#    p_cp = functions.skew_cp(x,mu,sigma,gamma)
#
##    print 'para injected            = ',e,w,a
##    print 'para with transformation = ',mu,sigma,gamma
##    print 'para recovered           = ',functions.transform_skew_para_from_cp_to_dp(mu,sigma,gamma)
#
#    mean_p = append(mean_p,stats.wmean(x,p))
#    p_corr = functions.skew(x,e-correct_mean_offset,w,a)
#    mean_p_corr = append(mean_p_corr,stats.wmean(x,p_corr))
##mod,solution,chi2 = functions.fit_skew_distribution(x,p,e,w)
#    mod,solution,chi2 = functions.fit_skew_distribution(x,p,e,w)
#    mean_fit_p = append(mean_fit_p,stats.wmean(x,mod))
#    #mean_fit_p2 = append(mean_fit_p2,solution[0][2] + functions.correct_mean_offset(solution[0][3],solution[0][4]))
#    mu_fitted = append(mu_fitted,solution[0])
#    gamma_fitted = append(gamma_fitted,solution[2])
#
#    subplot(211)
#    plot(x,p)
#    plot(x,p_cp,'r')
#    plot(x,mod,'g--')
#    ax = gca()
#    subplot(212,sharex = ax)
#    plot(x,p_corr)


xlabel('RV [m/s]')
ylabel('Flux')
savefig('figures/Profile_variation_with_alpha.pdf')

figure()
subplot(211)
title ('alpha varying, e fixed')
plot(skew_vect,mu_real,lw=12,c='y',label='mu real(analytical mean)')
plot(skew_vect,mean_p,lw=6,label='numerical mean')
#plot(skew_vect,mean_p_corr,lw=3,label='mean corrected')
plot(skew_vect,mean_fit_p,lw=2,c='r',label='numerical mean SN fitted')
plot(skew_vect,mu_fitted,lw=2,c='g',ls='--',label='fitted mean')
#plot(skew_vect,mean_skew_vect,lw=2,c='r',label='analytical mean')
xlabel(r'injected skewness para a')
ylabel('mean of the normal skew [m/s]')
legend(loc=1)


subplot(212)
plot(skew_vect,gamma_real,label='gamma real')
plot(skew_vect,gamma_fitted,lw=2,c='r',ls='--',label='gamma fitted')
#plot(skew_vect,mean_skew_vect,lw=2,c='r',label='analytical mean')
xlabel(r'injected skewness para a')
ylabel(r'fitted skewness para $\gamma$')
legend(loc=1)
savefig('figures/Mean_and_alpha_of_profile_variation.pdf')


figure()
plot(gamma_real,gamma_fitted-gamma_real)
ylabel(r'gamma_fitted-gamma_real')
xlabel(r'gamma_real')
savefig('figures/diff_between_real_and_fitted_gamma.pdf')





###########################################################################
############### Look at the chi2 close to 0
###########################################################################


#center,sigma,alpha = 10000.,5000.,0.0012
#y = functions.skew(x,center,sigma,alpha)
#y0 = 1./(sigma*sqrt(2*pi)) * exp(-(x-center)**2/(2.*sigma**2))
#figure()
#plot(x,y,lw=3)
#plot(x,y0,color='r',ls='--',lw=3)
#alpha_vect = arange(-0.1,0.1,0.0001)
#center_vect = arange(9950,10050,0.01)
#
#chi2 = array([sum((functions.skew(x,center0,sigma,alpha) - y)**2) for center0 in center_vect])
#figure()
#title('center = %.5f (real = %.5f)' % (center_vect[chi2.tolist().index(min(chi2))],center))
#plot(center_vect,chi2)
#
#chi2 = array([sum((functions.skew(x,center,sigma,alpha0) - y)**2) for alpha0 in alpha_vect])
#figure()
#title('alpha = %.5f (real = %.5f)' % (alpha_vect[chi2.tolist().index(min(chi2))],alpha))
#plot(alpha_vect,chi2)





#e = 0.0 # location
#w = 2.0 # scale
#
#x = linspace(-10,10+1e-10,n)
#
#
#fzz = skew(x,e,w,a) + norm.rvs(0,0.04,size=n) # fuzzy data
#
#def optm(l,x):
#    return skew(x,l[0],l[1],l[2]) - fzz
#
#best_sol = leastsq(optm,[0.1,0.1,0.1],(x,))
#
#print best_sol
#
#figure()
#title('best fit: %.2f,%.2f,%.2f' % (best_sol[0][0],best_sol[0][1],best_sol[0][2]))
#plot(x,fzz,'o')
#plot(x,skew(x,best_sol[0][0],best_sol[0][1],best_sol[0][2]),lw=3,c='r')
#plot(x,skew(x,best_sol[0][0],best_sol[0][1],0),lw=3,c='g')
#plot(x,norm.pdf(x,best_sol[0][0],best_sol[0][1]),lw=3,c='m',ls='--')
