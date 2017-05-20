
Hey there,

The LineProf code implements the line profile indicators described in Figueira 
et al. (2013, reference below) in python language, in a simple and easy-to-use 
manner, intuitive even for those non familiar with python, but allowing  
language users to expand on the provided algorithms. It was written by yours 
truly with valuable contributions from Jorge Humberto and João Faria, and tested
by willingfully compliant guinea pigs from CAUP. We will try to keep it updated 
on our (little) free time, so keep an eye on the site.

All the best,

P. (10 Feb 2014)
                                                    
DESCRIPTION
The LineProf code implements a series of line-profile analysis indicators and 
evaluates its correlation with RV data. The code receives as input a list of 
Cross-Correlation Functions and an optional list of associated RV. The code will 
evaluate the line-profile according to the indicators and compare it with the 
computed RV if no associated RV is provided, or with the provided RV otherwise.
The strength of the correlation is assessed by claculating the Pearson's 
correlation coefficient of the data and comparing it with that of an equivalent
uncorrelated dataset. Such dataset, or population of datasets is obtained by 
performing 100k times a Fisher-Yates shuffle of the indicator population. The 
correlation coefficient of the original dataset is compared with the shuffled 
one, and the z-score is provided. Assuming a gaussian distribution of the 
uncorrelated datasets correlation coefficients, one can associate to the z-score 
the one-sided probability of having such a value or larger drawn from an 
uncorrelated distribution. 
We considered as line profile indicators the BIS, including the different 
parametrizations BIS+ and BIS-, the indicators biGaussian and Vasy, as well as 
the FWHM. 

REQUIREMENTS
The code was extensively tested on python2.6 and 2.7, and even though it is 
expected to work properly on more recent versions, the user is advised to use 
one of these. The very standard libraries:

numpy
scipy
matplotlib    
pyfits

have to be installed and accecible. The program does not rely on any compiled 
code and the users are encourage to adapt the provided code to their one needs
and test it. It can also be used as a library in itself and incuded in larger 
code.

INPUT
The input to the program is the file "LineProfConfig.txt" and a list of the 
files to be processed, in ASCII format. The configuration parameters and file
properties are described below, using the enumeration 

  (1) - text file with list of observations to be processed. Each line 
	corresponds to the full path an observation to be processed. 

  (2) - type of files discriminated in observation list. The accepted types are
	"HARPS", "HARPS-N", "SOAP", "ASCII", and "rdb". The "HARPS" and 
	"HARPS-N"  are the standard output *_ccf_A_mask.fits as produced by the 
	DRS. The "SOAP" fits are produced with the SOAP code:

	http://adsabs.harvard.edu/abs/2012A%26A...545A.109B
	http://adsabs.harvard.edu/abs/2013A%26A...556A..19O

	The "ASCII" files are column and tab-separated files in which the first 
	column corresponds to the RV of each pixel and the second collumn the 
	corresponding CCF value. The "rdb" file is similar to "ASCII" with a 
	two-lines header, which are ignored in the reading process.

  (3) - flags used to toggle analysis of the different indicators. Values of 
	1/0 toggle the analysis on/off, respectively. The indicators are 
	decribed in detail in 

	http://adsabs.harvard.edu/abs/2013A%26A...557A..93F

  (4) - the RV to be compared with the Indicators is, by default, the center of 
	a Gaussian function fitted to the CCF.
	The user can overwrite this value by providing a FITS KW in which the RV 
	is stored in each file (as can be the case for HARPS-S/HARPS-N/SOAP 
	files) or an ASCII file with the RVs.
	The valid options, other than "None", are then the name of a 
	single-column .txt file with the RV in the same order as the observation 
	list or a KW full path (e.g. HIERARCH ESO DRS CCF RVC for HARPS).

A set of example parameters is already presented and and can be used along with 
the SOAP test data provided. A second HARPS data set is included for the user to 
explore. Take a look at "testdata/HARPS-HD224789/" and the observation list 
"testdata/HARPS-HD224789/HARPS_example.txt".
	
OUTPUT
The output is stored in the directory results and is composed of "Res" text 
files with the result of the indicators' analysis and correlation with RV's. A 
figure of Indicator vs. RV, with the best-fit slope is provided in the "plot" 
file, in .pdf format. A Wrap-up of all indicator's values is provided in the 
"WrapUp" text file.

TODO
Remove the annoying "Use the new widget gtk.Tooltip" warning present in some 
versions.

TROUBLESHOOTING
 - Did you make sure you are using python2.6/python2.7, and have installed and 
 accessible the libraries specified above?
