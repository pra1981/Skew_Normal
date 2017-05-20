import numpy as np
from scipy.optimize import leastsq

################################################################################
#
#		Definition of bigauss functions
#
###############################################################################

def F_bigauss(x, p):

	"""
	gaussian(x, A, mean, FWHM, cont)

	p[0] = A;
	p[1] = mean;
	p[2] = FWHM;
	p[3] = asym
	p[4] = cont.

	"""

	tau = -(4 * np.log(2) * (x - p[1])**2) / (np.abs(p[2]) * (1.0 + p[3]*np.sign(x-p[1])) )**2.0
	result = p[0] * np.exp( tau )+ p[4];
	return result;

def residuals_bigaussf(p, x, y, err=1.0):
	model = F_bigauss(x, p)
	return(y-model)/err

def bigaussfitter(x, y, p0, err=1.0):

	plsq = leastsq(residuals_bigaussf, p0, args=(x, y, err))
	#print "\n\tThe best-fit parameters for the bigaussian are: ", plsq
	return plsq[0]
