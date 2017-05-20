import string, sys
from astropy.io import fits
import numpy as np
from datetime import datetime
import time
import textwrap

import matplotlib.pyplot as plt
from matplotlib import rc
#set stuff for latex usage
rc('text', usetex=True)

plotting = True 	#set True/False to toggle on/off plotting capabilities"

################################################################################

def read_config(configfile):
	""" read LineProfile configuration file """
	# Get config parameters from file
	data = read_col(configfile)
	listname = [line[0] for line in data if line[-1] == '(1)'][0]
	filetype = [line[0] for line in data if line[-1] == '(2)'][0]
	selected_indicators = [float(line[0]) for line in data if line[-1] == '(3)']
	RVext = [string.split(string.join(line[:-1], " "), "#")[0] for line in data if line[-1] == '(4)'][0]

	# parameters check
	if (filetype not in ["HARPS", "HARPS-N", "SOAP", "ASCII", "rdb"]):
		print_line('WARNING', 'Type of files not recognized. Use one of the files described in the configuration file')
		sys.exit(0)

	outliers = [flag for flag in selected_indicators if flag not in [0,1]]
	if(outliers != []):
		print_line('WARNING', "Unrecognized selector value, use only 0/1.")
		sys.exit(0)


	print_line('Config', "Config file %s read successfully" % (configfile))
	print_line('Config', "Observations list is %s and filetype %s" % (listname, filetype))
	if(RVext=="None "):
		print_line('Config', "No external RV values were provided" )
	elif(RVext.endswith(".txt")):
		print_line('Config', "External RV file selected: %s" %(RVext))
	else:
		print_line('Config', "RV FITS Keyword selected: %s" %(RVext))

	return listname, filetype, selected_indicators, RVext

################################################################################

def read_obslist(listname, filetype, RVext):
	""" read dataset as written in listname of type filetype """

	print_line('Config', 'Reading data')

	#retrieve the files name from observations list
	list_of_files = read_col(listname)
	
	#creates a dummy vector to be replaced if RV values are provided 
	RVextvalues = []
	
	#Read the external RVs from a text file
	if(RVext.endswith(".txt")):
		RVextvalues  = read_floatcol(RVext)
	
	#reading the files and storing them in a list
	if(filetype in ["HARPS", "HARPS-N", "SOAP"]):
		#then fits will be used to read the files
		dataset = [read_obsFITS(obs_file[0], filetype) \
		for obs_file in list_of_files]
		if(RVext!="None " and RVext.endswith(".txt")==False):
			RVextvalues = [read_kw(obs_file[0], RVext) \
			for obs_file in list_of_files]	
	elif(filetype=="ASCII"):
		#two-collumn files will be used to read the files
		dataset = [read_2col(obs_file[0]) for obs_file in list_of_files]
	else:
		#then the file will have to be a rdb file.
		dataset = [read_2col_rdb(obs_file[0]) for obs_file in list_of_files]
	
	#Reading BJD if available		
	if(filetype=="HARPS"):
		#if filetype HARPS get the BJD as delivered by the pipeline
		BJDs = [read_kw(obs_file[0], "HIERARCH ESO DRS BJD") \
			for obs_file in list_of_files]
		print_line('Config', 'BJD red from Kw HIERARCH ESO DRS BJD')	
		print BJDs[0:4]
	elif(filetype=="HARPS-N"):
		#if filetype HARPS get the BJD as delivered by the pipeline
		BJDs = [read_kw(obs_file[0], "HIERARCH TNG DRS BJD") \
			for obs_file in list_of_files]
		print_line('Config', 'BJD red from Kw HIERARCH TNG DRS BJD')
	else:
		#No BJD available
		print_line('Config', 'BJD information not available')
		BJDs="None"	
		
	print_line('Config', 'Data reading complete')
	print "\n"
	
	#return list of file names, without path
	file_names = [obs_file[0].split("/")[-1] for obs_file in list_of_files]
	return [file_names, dataset, BJDs, RVextvalues]

################################################################################

def print_now(msg):
	""" Print to stdout (without carriage return) and flush right away.
	Useful to print in the same line """
	sys.stdout.write(msg)
	sys.stdout.flush()


################################################################################
def print_line(tag, msg):
	""" prints a formated line with a tag, message and time
	"""
	if len(msg) > 100:
		msg = textwrap.wrap(msg,width=100)
		string = '[{0:^10}] {1:.<100} [{2:^10}]\n'.format(tag, msg[0],datetime.now().strftime('%H:%M:%S'))
		for line in msg[1:]:
			string += '{0:^12} {1:.<96} {2:^12}\n'.format('', line ,'')
	else:
		string = '[{0:^10}] {1:.<100} [{2:^10}]\n'.format(tag, msg,datetime.now().strftime('%H:%M:%S'))

	sys.stdout.write(string)

################################################################################
#
#	Functions to read and write files in column-separated formats
#
#################################################################################
def read_col(filename):
	""" reads data from a file , ignoring lines started with the character
	"#". Returns a list with the blank-separated elements of each line """

	f = open(filename, "r")
	with open(filename, "r") as f:
		list_data = [line.split() for line in f.readlines() if not line.startswith('#')]

	return list_data

#################################################################################

def read_floatcol(filename):
	""" the same as read_col, but used to read one collumn of
	floating-point values and returning them """

	list_data = read_col(filename);

	col1 = [float(line[0]) for line in list_data]
	return [col1];

def read_2col(filename):
	""" the same as read_col, but used to read two collumns of
	floating-point values and returning them """

	list_data = read_col(filename);

	col1, col2 = [float(line[0]) for line in list_data],[float(line[1]) for line in list_data]

	return [col1, col2];

def read_2col_rdb(filename):
	""" the same as read_col, but used to read two columns of
	floating-point values and returning them, while ignoring the header
	(the first two lines) """

	list_data = read_col(filename);

	col1, col2 = [float(line[0]) for line in list_data[2:]],[float(line[1]) for line in list_data[2:]]

	return [col1, col2];

def read_7col(filename):
	""" the same as read_col, but used to read 7 string collumns """

	list_data = read_col(filename);

	col1, col2, col3, col4, col5, col6, col7 = [line[0] for line in list_data],\
	[line[1] for line in list_data], [line[2] for line in list_data],\
	[line[3] for line in list_data], [line[4] for line in list_data],\
	[line[5] for line in list_data], [line[6] for line in list_data]
	
	return [col1, col2, col3, col4, col5, col6, col7];

#################################################################################
def save_results(resultsdir, indTag, file_names, RV, vector2, Gauss_params,
	fitted_slope, slope_params, chi2, std_res, rPearson_data, dev_data, prob_data):
	"""Saves the data obtained by the analysis """

	#the time string is set once so that the two files have the same timestamp
	time_stamp = datetime.now().strftime('%Y-%m-%dT%Hh%Mm')
	FileName = '%(tag)s_%(time)s' %{'tag':indTag,'time':time_stamp}

	output_ASCII = '%(dir)sRes%(name)s.txt'%{'dir':resultsdir,'name':FileName}

	create_header(output_ASCII, indTag, slope_params, chi2, std_res, rPearson_data, dev_data, prob_data)
	append_2col_andpars(output_ASCII, file_names, RV, vector2, Gauss_params)

	#print "["+indTag+"] Results stored in ASCII format in "+output_ASCII+"."
	print_line(indTag, 'Results stored in ASCII format in %s' %output_ASCII)

	fig=plt.figure(1)
	plt.title(indTag+" Vs. RV")
	plt.xlabel("RV [km/s]")
	plt.ylabel(indTag+" [km/s]")
	plt.plot(RV, vector2, 'ro', label =r"$\rho$ = %.3f" % (rPearson_data))
	plt.plot(RV, fitted_slope, 'g-')
	plt.legend(loc='best')
	delta_x = np.max(RV) - np.min(RV)
	plt.xlim(np.min(RV)-0.1*delta_x, np.max(RV)+0.1*delta_x)
	
	delta_y = np.max(vector2) - np.min(vector2)
	plt.ylim(np.min(vector2)-0.1*delta_y, np.max(vector2)+0.1*delta_y)
	
	if(plotting):
		plt.show()

	fig.savefig('%(dir)splot%(name)s.pdf' %{'dir':resultsdir,'name':FileName}, facecolor='w', format='pdf',bbox_inches='tight')
	plt.close()

	print_line(indTag, 'Plot stored in %(file)s' %{'file':'%(dir)splot%(name)s.pdf' %{'dir':resultsdir,'name':FileName}})
	
	#return output_ASCII so that results can be merged later onto a single file
	return output_ASCII
	
#################################################################################
def create_header(output, indTag, slope_params, chi2, std_res, rPearson_data, dev_data, prob_data):
	""" creates header of output file"""
	with open(output, "w") as f:
		f.write("# Results of "+indTag+" analysis.\n")
		f.write("# The fitted slope was (m,b) (%f, %f) with a rchi2 of %f and a \
		res. st.d. of %f\n" % (slope_params[0], slope_params[1], chi2, std_res))
		f.write("# The RV and indicator dataset have an associated Pearson's r of %f, \
		being at %f sigma from an uncorrelated shuffled distribution.\n" \
			% (rPearson_data, dev_data))
		f.write("# The single-sided Gaussian probability of such an event is %f %%\n" % (prob_data*100.0))
		f.write("#\t filename \t\t RV \t "+indTag+"\t\tA(G)\tcenter(G)\tFWHM(G)\tcont(G)\n")

#################################################################################
def append_2col_andpars(output, file_names, data1, data2, Gauss_params):
	""" appends file_names, data1 and data2 vectors in a tab-separated ASCII file"""

	f = open(output, "a")
	data_to_file = zip(file_names, data1, data2, Gauss_params)

	[f.write("\t"+line[0]+"\t\t"+str(line[1])+"\t"+str(line[2])+"\t\t"+
	str(line[3][0])+"\t"+str(line[3][1])+"\t"+str(line[3][2])+"\t"+
	str(line[3][3])+"\n") for line in data_to_file]
	f.close();

def wrapRes(output_ASCII_list, BJDs, resultsdir):
	""" Collapses all results provided into one ASCII file in resultsdir"""
	
	#recovery of indtags from file
	indTag_list = [ (filename.split("Res")[-1]).split("_")[0] for filename in output_ASCII_list]
	
	time_stamp = datetime.now().strftime('%Y-%m-%dT%Hh%Mm')
	FileName = 'WrapUp_%(time)s' %{'time':time_stamp}

	output_WrapUp = '%(dir)s%(name)s.txt'%{'dir':resultsdir,'name':FileName}
		
		#reading the data
	Ind_data = read_7col(output_ASCII_list[0])[0:2], [read_7col(file_ind)[2] for file_ind in output_ASCII_list] 

	with open(output_WrapUp, "w") as f:
		f.write("# Wrap-up of analysis.\n")
		f.write("#\t filename \t\t")
		if(BJDs!="None"):
			#If BJDs are present add them after the name
			f.write(" BJD \t\t")
		f.write(" RV \t\t")	
		[f.write(indTag+"\t\t") for indTag in indTag_list]
		f.write("\n")
	
		for i in range(len(Ind_data[0][0])):
			#after some stackoverflow reading I was advised not to overuse comprehension lists
			f.write("\t"+Ind_data[0][0][i])	#writing filename
			if(BJDs!="None"):
				#If BJDs are present add them after the name
				f.write("\t\t"+("%.10f" % BJDs[i]))	#writing BJDs when available
			f.write("\t\t"+Ind_data[0][1][i])	#writing RV resulting from fitting
			for j in range(len(Ind_data[1])):
				f.write("\t\t"+Ind_data[1][j][i])	#writing available indicators
			f.write("\n")
	
	print '\n'
	print_line('WrapUp', 'Wrap-up of results stored in ASCII format in %s' %output_WrapUp)

#################################################################################
#
#	Functions to read FITS files
#
#################################################################################
def read_data(fitsfilename):
    """ read_data(fitsfilename):\n\
    #Read the FITS file and return an array containing the data. """

    filep = fits.open(fitsfilename,memmap=False);
    data = filep[0].data
    filep.close()
    return data

def read_kw(fitsfilename, kw):
    """ read the kw from fitsfilename and returns its value"""

    filep = fits.open(fitsfilename,memmap=False)
    head = filep[0].header[kw]
    filep.close()
    return head

#################################################################################
def read_obsFITS(obs_file, filetype):
    """ read the FITS file obs_file of the given type and return an
    RV table and CCF CCF
    """
    fits_image = fits.open(obs_file,memmap=False);
    data_st = fits_image[0].data
    head_st = fits_image[0].header
    fits_image.close()
    if(filetype=="HARPS"):
        #the sum CCF for HARPS is in the 73rd order
        sumCCF = data_st[72];
        #recover crval1 and cdelt1 to create vel table
        crval1 = head_st["CRVAL1"]
        cdelt1 = head_st["CDELT1"]
    elif(filetype=="HARPS-N"):
        #the sum CCF for HARPS-N is in the 70th order
        sumCCF = data_st[69];
        #overwritecrval1 and cdelt1
        crval1 = head_st["CRVAL1"]
        cdelt1 = head_st["CDELT1"]
    else:
        #if none of the options is selected, we have SOAP data
        sumCCF = data_st
        crval1 = -1*head_st["WIDTH0"]
        cdelt1 = head_st["STEP"]

    # velocity table creation
    vel_table = np.array([ crval1 + i*cdelt1 for i in range(len(sumCCF))])

    return [vel_table, sumCCF]

################################################################################
# Crazy stuff or what happens at 5 in the morning...
################################################################################

def print_head():			# (http://bigtext.org/)
	print 126*'#'
	print '{0:^120}'.format('\n'+'\
\t\t####### #     # #######    ### #     # ######  ###  #####     #    ####### ####### ######   #####  \n\
\t\t   #    #     # #           #  ##    # #     #  #  #     #   # #      #    #     # #     # #     # \n\
\t\t   #    #     # #           #  # #   # #     #  #  #        #   #     #    #     # #     # #       \n\
\t\t   #    ####### #####       #  #  #  # #     #  #  #       #     #    #    #     # ######   #####  \n\
\t\t   #    #     # #           #  #   # # #     #  #  #       #######    #    #     # #   #         # \n\
\t\t   #    #     # #           #  #    ## #     #  #  #     # #     #    #    #     # #    #  #     # \n\
\t\t   #    #     # #######    ### #     # ######  ###  #####  #     #    #    ####### #     #  #####  \n')
	print 126*'#'+'\n'

################################################################################

def print_foot():
	print '\n'
	print 126*'#'
	print '{0:^120}'.format('\n'+'\
	#######                     ###              #                     #######                             ### \n\
	   #    #    #   ##   ##### ###  ####       # #   #      #         #        ####  #      #    #  ####  ### \n\
	   #    #    #  #  #    #    #  #          #   #  #      #         #       #    # #      #   #  #      ### \n\
	   #    ###### #    #   #   #    ####     #     # #      #         #####   #    # #      ####    ####   #  \n\
	   #    #    # ######   #            #    ####### #      #         #       #    # #      #  #        #     \n\
	   #    #    # #    #   #       #    #    #     # #      #         #       #    # #      #   #  #    # ### \n\
	   #    #    # #    #   #        ####     #     # ###### ######    #        ####  ###### #    #  ####  ### \n')
	print 126*'#'+'\n'
