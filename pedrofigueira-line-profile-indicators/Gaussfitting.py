import numpy as np
from scipy.optimize import leastsq

import matplotlib.pyplot as plt
from matplotlib import rc
import time
import random as rnd
#set stuff for latex usage
rc('text', usetex=True)

from InOut import print_line

verbose = False	#set True/False to toggle talkative mode on/off
plotfit = False	#set True/False to toggle interactive plot of Gaussian fitting

#
#	Simple Gaussfitting function definitions
#

def F_gauss(x, p):

	"""
	gaussian(x, A, mean, FWHM, cont)

	p[0] = A;
	p[1] = mean;
	p[2] = FWHM;
	p[3] = cont.

	"""

	sigma = np.abs(p[2]) /( 2 * np.sqrt(2 * np.log(2)) );
	tau = -((x - p[1])**2) / (2*(sigma**2))
	result = p[0] * np.exp( tau )+ p[3];
	return result;

def residuals_gaussf(p, x, y, err=1.0):
	""" function that calculates residuals relative to a Gaussian fit"""
	model = F_gauss(x, p)
	return(y-model)/err

def example_fit():
	""" fit example """
	x = np.arange(100, dtype="double")
	p0 = [2.0, 2.2, 2.0, 2.5]
	err=[1.0]*100

	p_i = np.array([1.1,1.1,1.1,1.1])
	y = F_gauss(x, p_i);
	plsq = leastsq(residuals_gaussf, p0, args=(x, y, err))
	print "The best-fit parameters are: ", plsq

def gaussfitter(x, y, p0, err=1.0):

	plsq = leastsq(residuals_gaussf, p0, args=(x, y, err))
	
	if(verbose):
		print_line('Gauss', "The best-fit parameters are: %s" %plsq[0])
	
	#do modulus correction (note that modulus resulting from fitting
	#can be negative)
	results = [plsq[0][0], plsq[0][1], abs(plsq[0][2]), plsq[0][3]]
	return results
#
#	1D function definitions
#

def differentiate(y, x=None):
	""" function that differentiates a function using numpy diff """
	if(x==None):
		x=np.arange(0, len(y), 1.0);

	# Note that the differentiation reduces the vector dimension in 1
	# The other end is cut so that the result has the same center.

	numer=np.diff(y)
	deno=np.diff(x)

	result= numer[1:]/deno[1:];
	return result;

def peakfinder(y,x=None, smooth=0):
	""" function that finds the peak in a Gaussian function """

	if(x==None):
		x=np.arange(0, len(y), 1.0);

	# This peakfinder performs 2 derivatives and then extracts the maximum of the 2nd.
	# The first gets the variation of the spectra and finds local minima.
	deriv1=differentiate(y, x);
	if(smooth==1):
		# Smoothing in a size-3 box
		smoothed= np.zeros(len(deriv1))
		smoothed[0]=deriv1[0];
		smoothed[-1]=deriv1[-1];
		for i in range(1, len(deriv1)-1):
			smoothed[i]=np.mean(deriv1[i-1:i+2])
	else:
		smoothed=deriv1;

	# Peak listing
	zerolist=[];
	for i in range(len(smoothed)-1):
		if(smoothed[i]<0.0 and smoothed[i+1]>0.0):
			zerolist.append([x[i+2], y[i+1]]);
	if(zerolist==[]):
		print_line('Gauss', "No peaks were found")
		return zerolist;
	else:
		lower=zerolist[0][1];
		low_set=zerolist[0]
		# Here we assumed that the probability of finding 2 peaks with exactly the same value was negligible.
		for i in range(len(zerolist)):
			if(zerolist[i][1]<lower):
				lower=zerolist[i][1];
				low_set=zerolist[i]

		return low_set;


def fit_estimator(y,x=None, smooth=0):
	""" function that estimates starting parameters for the Gaussian fit """

	if(x==None):
		x=np.arange(0, len(y), 1.0);

	# Estimation starts by estimating the peak center
	peak = peakfinder(y,x, smooth);
	if (peak==[]):
		print_line('Gauss', "No peaks were found")
		return []

	center = peak[0]

	# Amplitude estimation by chopping the lower part of the sorted array
	sorted_y=np.sort(y);
	# 20 just looked good... no particular reason
	reduced_y=sorted_y[int(len(sorted_y)/20):]
	amplitude = (peak[1] - np.median(reduced_y));
	# Continuum is also estimated using the same principle
	continuum = np.median(reduced_y);

	# FWHM estimation by taking the points halfway between continuum and depth
	delta_x = x[1]-x[0] 	#we assume it is constant
	for i in range(0, len(x)):
		if (y[i] <= continuum + amplitude/2.0) :
			break

	FWHM=2.0*(peak[0]-x[i]);

	return [amplitude, center, FWHM, continuum];

def gauss_estimnfit(x, y, err=1.0):

	"""
	>>> import numpy as np
	>>> from gaussfit import *
	>>> x = np.arange(100, dtype="double")
	>>> p_i = np.array([-10.1,10.1,3.1,1.1])
	>>> y = F_gauss(x, p_i);
	>>> fit_estimator(y,x)
	[-7.1237720263844171, 10.0, 2.0, 1.1000000000000001]
	>>> gauss_estimnfit(x,y)
	"""

	p0 = fit_estimator(y,x);
	# Smoothing reaction?
	if(p0==[]):
		print "Trying to analyse with smoothing..."
		p_0 = fit_estimator(y,x, smooth=1);

	if(verbose):
		print "\t[Gauss] Estimation will be done using parameters:", p0;

	plsq = leastsq(residuals_gaussf, p0, args=(x, y, err))

	if(verbose):
		print "\n\t[Gauss]The best-fit parameters for the gaussian are: ", plsq

	#do modulus correction (note that modulus resulting from fitting
	#can be negative)
	results = [plsq[0][0], plsq[0][1], abs(plsq[0][2]), plsq[0][3]]
	return results

def dataset_gaussfit(dataset):
	""" Fits Gaussian function for each observation of the dataset """

	if(verbose==False):
		print_line('Gauss', 'Fitting Gaussian function to the data')

	Gauss_params = [gauss_estimnfit(obs[0], obs[1]) for obs in dataset]

	print_line('Gauss','Gaussian function fitting complete')

	if(plotfit):
		# plot fit of gaussian functions
		interGaussianplot(dataset, Gauss_params)

	return Gauss_params

#
#	Gaussian function iterative fitting
#

def interGaussianplot(dataset, Gauss_params, time_sleep=3):
	""" (mildly) interactive plotting """

	#old version
	#FittedData = [F_gauss(dataset[i][0], Gauss_params[i]) for i in range(len(dataset))]
	FittedData = [F_gauss(DataLine[0], GaussParam) for DataLine,GaussParam in zip(dataset,Gauss_params)]

	plt.figure(1)
	plt.ion()		# turns interactive mode on
	for DataLine, FittedDataLine in zip(dataset,FittedData):
		#Top plot
		plt.subplot(211)
		plt.ylabel("Counts on CCF [ ]")
		plt.title("CCF function with Gaussian fit.")
		plt.plot(DataLine[0], DataLine[1], 'ro', label="Obs")
		plt.plot(DataLine[0], FittedDataLine, 'g-', label ="Fit")
		plt.legend(loc='best')
		# Bottom Plot
		plt.subplot(212)
		plt.ylabel("OMC CCF [ ]")
		plt.xlabel("RV[km/s]")
		plt.plot(DataLine[0], DataLine[1]-FittedDataLine, 'k-')
		plt.draw()		# update window
		plt.show(False)
		time.sleep(time_sleep)	# in seconds
		plt.clf()
	plt.close()
	plt.ioff()		# turns interactive mode off
