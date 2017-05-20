import numpy as np
from scipy.optimize import leastsq
from random import randint
from scipy.stats import norm

from InOut import print_line

def r_Pearson(x, y):
	""" Calculus of the Pearson product-moment correlation coeficient
	http://en.wikipedia.org/wiki/Pearson_product-moment_correlation_coefficient
	"""

	x = np.array(x)
	y = np.array(y)

	x_mean = np.average(x)
	y_mean = np.average(y)

	r = np.sum( (x - x_mean) * (y - y_mean) ) / np.sqrt(np.sum((x - x_mean)**2) * np.sum((y - y_mean)**2))

	return r

def shuffle_FY(vector):
	""" Fisher-Yates shuffling according to Durstenfeldt/Knuth's algorithm
	http://en.wikipedia.org/wiki/Fisher-Yates """
	for index in xrange(len(vector)-1, 0, -1):
		ind_swap =randint(0, index)
		vector[ind_swap], vector[index] = vector[index], vector[ind_swap]
	return vector

def corr_MCcomp(indTag, vector1, vector2, MC_number=100000):
	"""calculates the significance of a correlation by shuffling the data
	to obtain non-correlated data and using it as H0 hypothesis"""
	r_origdata = r_Pearson(vector1, vector2)

	print_line(indTag,'Calculating significance of correlation using permutation test')

	vector_copy = list(vector2)
	# shuffle vector2 and calculate r_Pearson for each permutation
	r_shuffleddist = [r_Pearson(vector1, shuffle_FY(vector_copy)) for i in xrange(MC_number)] 
	#N.B.:The copy is shuffled repeatedly instead of creating a copy vector 
	# through the usage of list in order to gain time.
	print_line(indTag,'Permutation test complete')

	#This could be written in one line using the parameters
	# of the function CDF but is written in this way for clarity
	dev_data = abs(r_origdata)/np.std(r_shuffleddist)
	prob_data = (1.0-norm.cdf(dev_data))

	print_line(indTag,"The Pearson's rho is %(var1)f, and is at %(var2).3f sigma and thus \
having a prob of %(var3).3f %% when assuming a (single-tailed) Gaussian distribution." \
			%{'var1':r_origdata, 'var2':dev_data,'var3':(prob_data*100.0) })

	#print "["+indTag+"] The Pearson's rho is %(var1)f, and is at %(var2).3f sigma and thus \
	#having a prob of %(var3).3f %% when assuming a (single-tailed) Gaussian distribution." \
	#%{'var1':r_origdata, 'var2':dev_data,'var3':(prob_data*100.0)}

	return [r_origdata, dev_data, prob_data]

################################################################################
#
#		Linear fit functions
#
###############################################################################

def polyn1dg(x, p):
        """ fitting a slope """
        result = p[0]*x + p[1];
        return result;

def residuals_pol1(p, x, y, err=1.0):
	""" residuals after fitting a slope """
        model = polyn1dg(x, p)
	return (y-model)/err

def fitter_pol1(p0, x, y, err):
	""" function that fits slopes to data """
	plsq = leastsq(residuals_pol1, p0, args=(x, y, err))
	#print "The best-fit parameters are: ", plsq
	return plsq[0]

def fitter_slope(x,y):
	"""Functions that fits a slope assuming as a starting point a negative
	1:1 correlation (LM will do )"""
	#slope fitting part
	p0 = [-1.0, 1.0];
	err = [1.0] * len(x)
	p0f = fitter_pol1(p0, x, y, err);
	fitted_line = polyn1dg(x, p0f)

	#calculating the chi2 and the std of the residuals
	# http://en.wikipedia.org/wiki/Standard_deviation
	residuals2 = ((y - fitted_line))**2
	std_res = np.sqrt((1.0/(len(residuals2)-1.0))*np.sum(residuals2));

	#the reduced chi^2 is calculated using the Person's chi squared test
	# http://en.wikipedia.org/wiki/Goodness_of_fit
	rchi2 = 1.0/(len(y)-len(p0f)) * np.sum(residuals2)

	#print("The final red chi^2 for the solution is %12.5f." % rchi2)

	return [p0f, fitted_line, rchi2, std_res]

### small test function
def _test_fitter_slope():
	x = np.array([1,2,3,4,5,6])
	y = np.array([12, 24, 36, 48, 60, 75])
	p, fitted, rchi2, std_res = fitter_slope(x, y)
	assert np.allclose(p, [12.42857143, -1.])
	assert np.allclose([rchi2, std_res], [1.0714285714285714, 0.92582009977255153])
	# test passed!
