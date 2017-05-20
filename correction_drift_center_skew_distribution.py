#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib.pyplot as pl


from pylab import *
from scipy.special import erf
from scipy.stats import norm
from scipy.optimize import leastsq

from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm
import itertools

import sys
sys.path.append('/Users/xav/uni/phd/program/functions/')

import stats


########################################################
#### Functions
########################################################

# For more information about the skew distribution, see "http://azzalini.stat.unipd.it/SN/Intro/intro.html"
#and have a look at "test/skew_distribution_for_CCF"
def pdf(x):
    return 1/sqrt(2*pi) * exp(-x**2/2)

def cdf(x):
    return (1 + erf(x/sqrt(2))) / 2

def skew(x,e=0,w=1,a=0):
    t = (x-e) / (1.*w)
    return 2./w * pdf(t) * cdf(a*t)
## You can of course use the scipy.stats.norm versions
## return 2 * norm.pdf(t) * norm.cdf(a*t)


def calculate_mean_normal_kew(a,w):
    
    mean_skew = e + w * sqrt(2./pi) * a/(sqrt(1+a**2))
    
    ampli = 2e3
    p = ampli * skew(x,e-mean_skew,w,a)
    return stats.wmean(x,1-p)


########################################################
#### Main
########################################################

n = 2**10

x = linspace(-50000,50000+1e-10,n)

#vectors for skew and sigma
skew_vect = arange(-1,1+1e-10,0.01)
sigma_step = 1000.
sigma_vect = arange(sigma_step,5000+sigma_step/2.,sigma_step)

#Create the mesh grid and calculate the value of the mode for each point
A, W = np.meshgrid(skew_vect,sigma_vect)
zs = array([calculate_mean_normal_kew(a,w) for a,w in zip(np.ravel(A), np.ravel(W))])
M = zs.reshape(A.shape)


fig = figure()
ax = fig.gca(projection='3d')
ax.plot_surface(A, W, M, rstride=8, cstride=8, alpha=0.3)
cset = ax.contourf(A, W, M, zdir='z', offset=-200, cmap=cm.coolwarm)
cset = ax.contourf(A, W, M, zdir='x', offset=-1.5, cmap=cm.coolwarm)
cset = ax.contourf(A, W, M, zdir='y', offset=-1, cmap=cm.coolwarm)

ax.set_xlabel('A')
#ax.set_xlim(-40, 40)
ax.set_ylabel('W')
#ax.set_ylim(-40, 40)
ax.set_zlabel('Mean')
#ax.set_zlim(-100, 100)

savefig('figures/3d_plot_mean_corrected.pdf')
