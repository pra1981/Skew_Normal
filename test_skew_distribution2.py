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

e = -2.33787593e+04 # location
w = 2.74078830e+03 # scale

x = linspace(-50000,50000+1e-10,n)

skew_vect = arange(-0.3,0.3+1e-10,0.002)
mean_p,mean_p_corr,mean_fit_p,mean_skew_vect,e_fitted,a_fitted = [],[],[],[],[],[]
for i,a in enumerate(skew_vect):
    
    correct_mean_offset = functions.correct_mean_offset(w,a)
    
    if i%1000 == 0:
        print i,' over ',len(skew_vect)
    cte = 1.23212844e+08
    ampli = 3.88982676e+11
    p = cte-ampli*functions.skew(x,e,w,a)
    mean_p = append(mean_p,stats.wmean(x,cte-p))
    p_corr = cte-ampli*functions.skew(x,e-correct_mean_offset,w,a)
    mean_p_corr = append(mean_p_corr,stats.wmean(x,cte-p_corr))
    mean_skew_vect = append(mean_skew_vect,e+correct_mean_offset)
    mod,solution,init_mod,chi2,chi2_red = functions.fit_skew_distribution_for_CCF(x,p,3000*cte,e,2*sqrt(2.*log(2.))*w)
    mean_fit_p = append(mean_fit_p,stats.wmean(x,solution[0]-mod))
    #mean_fit_p2 = append(mean_fit_p2,solution[0][2] + functions.correct_mean_offset(solution[0][3],solution[0][4]))
    e_fitted = append(e_fitted,solution[2])
    a_fitted = append(a_fitted,solution[4])
    
    subplot(211)
    plot(x,p)
    plot(x,mod)
    subplot(212)
    plot(x,p_corr)


#skew_vect = arange(0,0.2+1e-10,0.001)
#mean_p,mean_p_corr,mean_fit_p,mean_skew_vect,e_fitted,a_fitted = [],[],[],[],[],[]
#for i,a in enumerate(skew_vect):
#    
#    correct_mean_offset = functions.correct_mean_offset(w,a)
#    
#    if i%1000 == 0:
#        print i,' over ',len(skew_vect)
#    p = functions.skew(x,e,w,a)
#    mean_p = append(mean_p,stats.wmean(x,p))
#    p_corr = functions.skew(x,e-correct_mean_offset,w,a)
#    mean_p_corr = append(mean_p_corr,stats.wmean(x,p_corr))
#    mean_skew_vect = append(mean_skew_vect,e+correct_mean_offset)
#    mod,solution = functions.fit_skew_distribution(x,p,e,w)
#    mean_fit_p = append(mean_fit_p,stats.wmean(x,mod))
#    #mean_fit_p2 = append(mean_fit_p2,solution[0][2] + functions.correct_mean_offset(solution[0][3],solution[0][4]))
#    e_fitted = append(e_fitted,solution[0])
#    a_fitted = append(a_fitted,solution[2])
#    
#    subplot(211)
#    plot(x,p)
#    plot(x,mod)
#    subplot(212)
#    plot(x,p_corr)






xlabel('RV [m/s]')
ylabel('Flux')
savefig('figures/Profile_variation_with_alpha.pdf')

figure()
plot(skew_vect,mean_skew_vect,lw=8,c='y',label='analytical mean distri')
plot(skew_vect,mean_p,lw=5,label='mean distri')
plot(skew_vect,mean_p_corr,lw=5,label='mean distri corr')
plot(skew_vect,mean_fit_p,lw=2,c='r',label='mean distri fit')
plot(skew_vect,e_fitted,lw=2,c='r',ls='--',label='e fitted')
#plot(skew_vect,mean_skew_vect,lw=2,c='r',label='analytical mean')
xlabel(r'skewness para $\alpha$')
ylabel('Mean of the profile [m/s]')
legend()
savefig('figures/Mean_of_profile_variation_with_alpha.pdf')

figure()
plot(skew_vect,(e_fitted-e)/e)
plot(skew_vect,a_fitted/mean(skew_vect),lw=3)




###########################################################################
############### Look at the chi2 close to 0
###########################################################################


center,sigma,alpha = 10000.,5000.,0.01
y = functions.skew(x,center,sigma,alpha)
y0 = 1./(sigma*sqrt(2*pi)) * exp(-(x-center)**2/(2.*sigma**2))
figure()
plot(x,y,lw=3)
plot(x,y0,color='r',ls='--',lw=3)
alpha_vect = arange(-0.1,0.1,0.0001)
center_vect = arange(9950,10050,1.)

chi2 = array([sum((functions.skew(x,center0,sigma,alpha) - y)**2) for center0 in center_vect])
figure()
title('center = %.5f (real = %.5f)' % (center_vect[chi2.tolist().index(min(chi2))],center))
plot(center_vect,chi2)

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
