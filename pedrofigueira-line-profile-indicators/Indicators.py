import numpy as np
import sys

import matplotlib.pyplot as plt
from matplotlib import rc

#set stuff for latex usage
rc('text', usetex=True)

import StatTools
from InOut import save_results, print_line
import Gaussfitting, BiGaussfitting

################################################################################
def run_indicator_analysis(indTag, dataset, Gauss_params, file_names, resultsdir, RVextvalues):
	""" do BIS analysis and store output in resultsdir """
	print '\n'
	print_line(indTag, 'Starting %(Tag)s analysis' %{'Tag':indTag})

	#Definition of the indicator to use from the indTag parameter provided
 	IndicatorDict = {"BIS":bissector,
 	                 "BIS-":bissector,
 	                 "BIS+":bissector,
 	                 "biGauss":biGauss,
 	                 "Vasy":Vasy, 
 	                 "Vspan":Vspan,
 	                 "FWHM": FWHM_analysis}
	indicator = IndicatorDict[indTag]

	#association of flux level to the different BIS 
 	FluxLevelsDict = {"BIS":[0.1, 0.4, 0.6, 0.9],"BIS-":[0.3, 0.4, 0.6, 0.7],"BIS+":[0.1, 0.2, 0.8, 0.9],"biGauss":None,"Vasy":None, "Vspan":None}

	if(indicator == bissector):
		#use flux levels in case bisector is selected
		Ind_pars = [indicator(DataLine[0],DataLine[1], GaussParam, FluxLevelsDict[indTag]) for DataLine,GaussParam in zip(dataset,Gauss_params)]
	else:
		Ind_pars = [indicator(DataLine[0],DataLine[1], GaussParam) for DataLine,GaussParam in zip(dataset,Gauss_params)]

	print_line(indTag, '%(Tag)s analysis complete' %{'Tag':indTag})

	RV, Ind = rearrange_indicator_data(Ind_pars)
	
	if(RVextvalues!=[]):
		#RV values to be used become RV external values
		print_line(indTag, 'Note: External RV values will be used.')
		RV = np.array(RVextvalues)
		
	slope_params, fitted_slope, chi2, std_res = StatTools.fitter_slope(RV, Ind);
	rPearson_data, dev_data, prob_data = StatTools.corr_MCcomp(indTag, RV, Ind)

	output_ASCII = save_results(resultsdir, indTag, file_names, RV, Ind, 
		                        Gauss_params, fitted_slope, slope_params, 
		                        chi2, std_res, rPearson_data, dev_data, prob_data)
	
	return output_ASCII
	
################################################################################
def bissector(RV, CCF, params, FluxLevels, plot = 0):
	# If plot=1: plot bissector as a function of depth

	FluxTopMin,FluxTopMax,FluxBottomMin,FluxBottomMax = FluxLevels

	depth, bis_b, bis_r, v0 = CCF_interpolate(RV, CCF, params)
	bis = (bis_b+bis_r)/2.

	qq = np.greater_equal(depth,FluxTopMin)*np.less_equal(depth,FluxTopMax)
	RV_top = np.average(np.compress(qq,bis))
	qq = np.greater_equal(depth,FluxBottomMin)*np.less_equal(depth,FluxBottomMax)
	RV_bottom = np.average(np.compress(qq,bis))
	span = RV_top-RV_bottom

	if(plot == 1):
		plt.figure(1)
		plt.plot(bis, (1.0-depth), 'ro')
		plt.title("BIS")
		plt.xlabel("BIS [km/s]")
		plt.ylabel("1-depth [ ]")
		plt.show()
		plt.close()
	return depth,bis,span,v0

#################################################################################
def biGauss(RV, CCF, Gauss_params, plot = 0 ):
	""" Function that calculates the BiGaussian value"""
	#LM fitting for a bigaussian taking as departing point the gaussian and assymetry=0
	dep_params = list(Gauss_params)
	dep_params = np.array([dep_params[0], dep_params[1], dep_params[2], 0.0, dep_params[3]])
	biG_params = BiGaussfitting.bigaussfitter(RV, CCF, dep_params)

	fitted =  BiGaussfitting.F_bigauss(RV, biG_params);

	if(plot == 1):
		plt.figure(1)
		#Top plot
		plt.subplot(211)
		plt.title("CCF function with bigaussian fit")
		plt.ylabel("Counts on CCF [ ]")
		plt.plot(RV, CCF, 'ro', label="Obs")
		plt.plot(RV, fitted, 'g-', label ="Fit")
		plt.legend(loc='best')
		#Bottom plot
		plt.subplot(212)
		plt.ylabel("OMC CCF [ ]")
		plt.xlabel("RV[km/s]")
		plt.plot(RV, CCF-fitted, 'k-')
		plt.show()
		plt.close()

	A_biG, center_biG, FWHM_biG, asym_biG, c_biG = biG_params
	center_G = Gauss_params[1]
	DeltaRV = center_biG - center_G

	return [center_biG, asym_biG, DeltaRV, center_G]

def Vasy(RV, CCF, Gauss_params, plot = 0):
	""" Function that calculates the vasy asymetry indicator """
	center_G = Gauss_params[1]
	dummy1, dummy2 = 0, 0	# created as dummys to obtain the same number of outputs

	depth, lambda_b, lambda_r, v0 = CCF_interpolate(RV, CCF, Gauss_params)

	# calculate weight of the left wing and the right one
	weight_r = weight_bouchy(lambda_r, depth)
	weight_b = weight_bouchy(lambda_b, depth)

	average_weight = (weight_r + weight_b) / 2.0
	
	Vasy_res = np.sum((weight_r - weight_b) * average_weight) / np.sum(average_weight)
	
	return [dummy1, dummy2, Vasy_res, center_G]

def Vspan(RV, CCF, Gauss_params, plot = 0):
	""" Function that calculates the Vspan indicator"""
	center_G = Gauss_params[1]
	dummy1, dummy2 = 0, 0	# created as dummys to obtain the same number of outputs

	#LM fitting for a gaussian but considering only the top
	left_limit = Gauss_params[1] - Gauss_params[2] /(2 * np.sqrt(2 * np.log(2)))
	right_limit = Gauss_params[1] + Gauss_params[2] /(2 * np.sqrt(2 * np.log(2)))
	
	Cindex = np.less(RV, left_limit) + np.greater(RV, right_limit)
	RV_sel = np.compress(Cindex, RV)
	CCF_sel = np.compress(Cindex, CCF)
	dep_params = list(Gauss_params)
	
	ccf_params_top = Gaussfitting.gaussfitter(RV_sel, CCF_sel, dep_params)

	#LM fitting for a gaussian but considering only the bottom
	left_limit_3s = Gauss_params[1] - 3.0*Gauss_params[2] /(2 * np.sqrt(2 * np.log(2)))
	right_limit_3s = Gauss_params[1] + 3.0*Gauss_params[2] /(2 * np.sqrt(2 * np.log(2)))
	
	left_limit = Gauss_params[1] - Gauss_params[2] /(2 * np.sqrt(2 * np.log(2)))
	right_limit = Gauss_params[1] + Gauss_params[2] /(2 * np.sqrt(2 * np.log(2)))

	Cindex = np.less(RV, left_limit_3s) + np.greater(RV, left_limit)*np.less(RV, right_limit) + np.greater(RV, right_limit_3s)
	
	RV_sel = np.compress(Cindex, RV)
	CCF_sel = np.compress(Cindex, CCF)
	dep_params = list(Gauss_params)
	
	ccf_params_bottom = Gaussfitting.gaussfitter(RV_sel, CCF_sel, dep_params)
	
	#Vspan is the Rv difference between teh top and the bottom
	Vspan_res = ccf_params_top[1] - ccf_params_bottom[1]
	
	return [dummy1, dummy2, Vspan_res, center_G]

def FWHM_analysis(RV, CCF, Gauss_params, plot = 0):
	""" Function that calculates the vasy asymetry indicator """
	center_G = Gauss_params[1]
	FWHM = Gauss_params[2]
	dummy1, dummy2 = 0, 0	# created as dummys to obtain the same number of outputs
	
	return [dummy1, dummy2, FWHM, center_G]

################################################################################

def rearrange_indicator_data(Ind_pars):
	""" Function that re-arranges data and returns numpy arrays """

	RV = np.array([params[3] for params in Ind_pars])
	Ind = np.array([params[2] for params in Ind_pars])

	return [RV, Ind]

################################################################################

def CCF_interpolate(RV, CCF, params):
	# Bissector calculus using Chris-like method of squared interpolation
	# The calculus itself is based on Chris variables and code

	k, v0, fwhm, c = params

	sigma = fwhm/2./np.sqrt(2.*np.log(2.))
	norm_CCF = -c/k*(1.-CCF/c)
	nstep = 100
	margin = 5
	depth = np.arange(nstep-2*margin+1,dtype=float)/nstep + float(margin)/nstep

	# mean RV for each segment of the CCF
	MeanRV = [(RV[i]+RV[i+1])/2. for i in range(len(CCF)-1)]

	# derivatives for each segment of the CCF
	ExpV= [np.exp(-(v-v0)**2/2/sigma**2)/sigma**2 for v in MeanRV]
	dCCFdRV = [-(v-v0)*expV for v,expV in zip(MeanRV,ExpV)]
	d2CCFdRV2 = [((v-v0)**2/sigma**2-1)*expV for v,expV in zip(MeanRV,ExpV)]
	d2RVdCCF2 = np.array([-d1/d2**3 for d1,d2 in zip(d2CCFdRV2,dCCFdRV)])
	
	# not-null range (a ver como simplificar)
	iRange = [ii for ii in range(len(CCF)-1) if (max(norm_CCF[ii],norm_CCF[ii+1]) >= depth[0]) & (min(norm_CCF[ii],norm_CCF[ii+1]) <= depth[-1])]

 	# parameters ?? for each segment of the CCF
 	p = np.zeros([len(CCF),3],'d')
 	p[iRange,2] = np.array(d2RVdCCF2[iRange])/2.
 	p[iRange,1] = np.array([(RV[i+1]-RV[i]-p[i,2]*(norm_CCF[i+1]**2-norm_CCF[i]**2))/(norm_CCF[i+1]-norm_CCF[i]) for i in iRange])
 	p[iRange,0] = np.array([RV[i]-p[i,1]*norm_CCF[i]-p[i,2]*norm_CCF[i]**2 for i in iRange])

	# Indexes where "norm_CCF > dd"
	Indexes = np.array([[np.where(norm_CCF > dd)[0][0]-1,np.where(norm_CCF > dd)[0][-1]] for dd in depth])
	IndexBlue, IndexRed = [ind[0] for ind in Indexes],[ind[-1] for ind in Indexes]

 	# Bisector definition
 	bis_b = np.array([p[i_b,0]+p[i_b,1]*dd+p[i_b,2]*dd**2 for i_b,dd in zip(IndexBlue,depth)])
 	bis_r = np.array([p[i_r,0]+p[i_r,1]*dd+p[i_r,2]*dd**2 for i_r,dd in zip(IndexRed,depth)])
	
	"""
	#Compare
	old_p = np.zeros([len(CCF),3],'d')
	old_bis_b = np.zeros(len(depth),'d')
	old_bis_r = np.zeros(len(depth),'d')
	for i in range(len(CCF)-1):
		if (max(norm_CCF[i],norm_CCF[i+1]) >= depth[0]) & (min(norm_CCF[i],norm_CCF[i+1]) <= depth[-1]):
			v = (RV[i]+RV[i+1])/2.
			old_dCCFdRV = -(v-v0)/sigma**2*np.exp(-(v-v0)**2/2/sigma**2)
			old_d2CCFdRV2 = ((v-v0)**2/sigma**2-1)/sigma**2*np.exp(-(v-v0)**2/2/sigma**2)
			old_d2RVdCCF2 = -old_d2CCFdRV2/old_dCCFdRV**3
			#print "ders:", old_dCCFdRV, dCCFdRV[i], old_d2CCFdRV2, d2CCFdRV2[i], old_d2RVdCCF2, d2RVdCCF2[i]
			old_p[i,2] = old_d2RVdCCF2/2
			old_p[i,1] = (RV[i+1]-RV[i]-old_p[i,2]*(norm_CCF[i+1]**2-norm_CCF[i]**2))/(norm_CCF[i+1]-norm_CCF[i])
			old_p[i,0] = RV[i]-old_p[i,1]*norm_CCF[i]-old_p[i,2]*norm_CCF[i]**2
			#print "p's[2]:", old_p[i,2],p[i,2], "\np's[1]", old_p[i,1], p[i,1], "\np's[0]",old_p[i,0],p[i,0]
	for j in range(len(depth)):
		i_b = norm_CCF.argmax()
		while (norm_CCF[i_b] > depth[j]) & (i_b > 1): i_b = i_b-1
		i_r = norm_CCF.argmax()
		while (norm_CCF[i_r+1] > depth[j]) & (i_r < len(CCF)-2): i_r = i_r+1
		print "\n j for depth", j, "\ni_b", i_b, IndexBlue[j], "\ni_r", i_r, IndexRed[j]
		old_bis_b[j] = p[i_b,0]+p[i_b,1]*depth[j]+p[i_b,2]*depth[j]**2
		old_bis_r[j] = p[i_r,0]+p[i_r,1]*depth[j]+p[i_r,2]*depth[j]**2
		print "bis_b:", old_bis_b[j], bis_b[j], "\nbis_r", old_bis_r[j], bis_r[j]
	"""
	
	return depth, bis_b, bis_r, v0

def weight_bouchy(lambdas, flux):
	#calculates the weight as a function of the pixel using the formula 8 of Bouchy, Pepe & Queloz (2001)
	dA_over_dlambda = (flux[1:] - flux[:-1]) / (lambdas[1:] - lambdas[:-1])
	#add one value in the end
	dA_over_dlambda = np.append(dA_over_dlambda, 0.0)
	
	weight = lambdas*lambdas*dA_over_dlambda*dA_over_dlambda/flux
	
	return weight

################################################################################
def run_indicator_analysis_CUSTOM(indTag, dataset, Gauss_params, file_names, resultsdir):
	""" do BIS analysis and store output in resultsdir 
	this analysis covers the hypothesys of having BIScustom, an advanced 
	feature capability"""
	
	print '\n'
	print_line(indTag, 'Starting %(Tag)s analysis' %{'Tag':indTag})

 	IndicatorDict = {"BIS":bissectorAll,"BIS-":bissectorAll,"BIS+":bissectorAll,"BIScustom":bissectorAll,"biGauss":biGauss,"Vasy":Vasy}
	indicator = IndicatorDict[indTag]


 	FluxLevelsDict = {"BIS":[0.1, 0.4, 0.6, 0.9],"BIS-":[0.3, 0.4, 0.6, 0.7],"BIS+":[0.1, 0.2, 0.8, 0.9],"biGauss":None,"Vasy":None}
 	if (indTag=="BIScustom"):
 		FluxBottomMin,FluxBottomMax,FluxTopMin,FluxTopMax = 1,1,1,1
 		while not FluxBottomMin < FluxBottomMax < FluxTopMin < FluxTopMax:
 			print_line(indTag,'Please enter custom flux limits (FluxBottomMin < FluxBottomMax < FluxTopMin < FluxTopMax)')
 			try:
				FluxBottomMin,FluxBottomMax = input('{:^12} {:.<19}'.format('','VBottom (min,max): '))
				FluxTopMin,FluxTopMax = input('{:^12} {:<16}'.format('','Vtop (min,max): '))
				FluxLevelsDict["BIScustom"] = [FluxBottomMin,FluxBottomMax,FluxTopMin,FluxTopMax]
			except GeneratorExit,SystemExit:
				sys.exit()
			except:
				print 'Wrong input!'

	try:
		Ind_pars = [indicator(DataLine[0],DataLine[1], GaussParam, FluxLevelsDict[indTag]) for DataLine,GaussParam in zip(dataset,Gauss_params)]
	except UnboundLocalError:
		Ind_pars = [indicator(DataLine[0],DataLine[1], GaussParam) for DataLine,GaussParam in zip(dataset,Gauss_params)]

	print_line(indTag, '%(Tag)s analysis complete' %{'Tag':indTag})

	RV, Ind = rearrange_indicator_data(Ind_pars)

	slope_params, fitted_slope, chi2, std_res = StatTools.fitter_slope(RV, Ind);
	rPearson_data, dev_data, prob_data = StatTools.corr_MCcomp(indTag, RV, Ind)

	save_results(resultsdir, indTag, file_names, RV, Ind, Gauss_params, fitted_slope, slope_params,
		chi2, std_res, rPearson_data, dev_data, prob_data)

################################################################################	
	
