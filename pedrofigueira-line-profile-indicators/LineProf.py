# LineProf Code by Pedro Figueira
# with contributions from Jorge Humberto & Joao Faria.
# Last update February 2014, keep an eye on:
# https://bitbucket.org/pedrofigueira/line-profile-indicators

import numpy as np
import os
import textwrap
try:
	os.system('clear')
except:
	os.system('cls')
	
# user-defined libraries
import InOut, Indicators, Gaussfitting

configfile = "LineProfconfig.txt"		# LineProf configuration file
resultsdir = "results/"				# results destination directory

#################################################################################
def run_analysis(configfile, resultsdir):
	""" master function that runs line profile analysis of the dataset """
	
	# Header
	InOut.print_head()

	#reading config from configfile
	listname, filetype, selected_indicators, RVext = InOut.read_config(configfile)

	#read observations from obslist
	file_names, dataset, BJDs, RVextvalues = InOut.read_obslist(listname, filetype, RVext)

	#fit Gaussian function to the data
	Gauss_params = Gaussfitting.dataset_gaussfit(dataset)

	# Identify selected indicators
	IndicatorList = ["BIS","BIS-","BIS+","biGauss","Vasy", "Vspan", "FWHM"]
	IndSelected = [ind for ind,selected in zip(IndicatorList, selected_indicators) if selected == 1]
	
	#Apply indicators
	output_ASCII_list = [Indicators.run_indicator_analysis(ind, dataset, Gauss_params, file_names, resultsdir, RVextvalues) for ind in IndSelected]
	
	#Wrap the results in one single file
	InOut.wrapRes(output_ASCII_list, BJDs, resultsdir)	
	
	# The End
	InOut.print_foot()

#################################################################################

if __name__ == "__main__":
	# testing function working
	run_analysis(configfile, resultsdir)


