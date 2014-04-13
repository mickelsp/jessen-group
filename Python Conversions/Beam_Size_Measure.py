## Overview of Beam_Size_Measure
#This program takes the results of knife edge measurements of laser beams, which should have the form of 
#an error function if the beam is Gaussian, and turns them into a 1/e^2 beam size.
#The shortcut is to take only the 90% of max power and the 10% of max power positions.

import numpy as np # Import numpy functions
import matplotlib.pyplot as plt # Import matplotlib functions
import scipy.special as sp # Import scipy functions
import scipy.optimize as optimize # Import optimize function for least squares fitting

plt.close()

## Options
fontsize = 16
shortcut1090 = 1 # 1 if you want to use positions of only 90% and 10% maximum power to determine beam size (shortcut for beam profiling)

if shortcut1090 == 1:
	## Import data for 90/10 measurements
	directoryname = '/Users/work/Documents/Arizona/Nanofibers/Laser Profiling/Probe Beam Profiling/Before Fiber Link/'
	inputfilename = 'ProbeLaser_FiberLink_ReverseOutput_20131008.dat'
	outputfilename = 'ProbeLaser_FiberLink_ReverseOutput_20131008.txt'
	filename = directoryname + inputfilename
	
	measurementposition = [] 
	position90=[]
	position10=[] # Initialize vectors
	inputfile = open(filename, 'r') # Open file for reading
	for line in inputfile:
		line = line.strip() 			# Read in a line at a time
		if not line.startswith("%"): 	# Only read a line if it doesn't start with "%", the comment symbol used in Matlab.
			columns = line.split()
			measurementposition.extend([float(columns[0])]) # Append value of measurementposition from this line to the list of measurement position values
			position90.extend([float(columns[1])])
			position10.extend([float(columns[2])])
	inputfile.close()

	## Convert lists into numeric vectors
	measurementposition = np.array(measurementposition) #[cm]
	position90 = np.array(position90) #[mm]
	position10 = np.array(position10) #[mm]

	## Calculate waist size based on Siegman formula (see p. 94 Pascal's Lab Notebook 1)
	waistsize = abs(position90-position10)/1.28 #[mm] 1/e^2 beam radius as a function of measurement position along length of beam
	
	## Plot waist size versus distance from lens]
	plt.figure(1)
	plt.plot(measurementposition,waistsize,'ro')
	plt.show()
	plt.xlabel('Longitudinal Position [cm]')
	plt.ylabel('Waist Size [mm]')
	#xlim([0 1.13])
	#ylim([0 max(waistsize(:))+0.1.*max(waistsize(:))])
	#legend('Vertical Size','Horizontal Size','Location','best')
	plt.title('Beam Sizes of 20 Hz NdYAG Laser')
	
elif shortcut1090 == 0:
	## Analyze full knife edge data set
	measurementposition = [34.5] #[cm] Position relative to beam's zero point
	directoryname = 'Z:/Corrosion Removal/Phase I/Data/2014.04.04 NdYAG 20 Hz Beam Profile/'
	filenames = ['NdYAG20Hz_BeamProfile_20140404_Full.txt','NdYAG20Hz_BeamProfile_20140404_Full.txt']
	numfiles = len(measurementposition) # There should be one data file per measurement position

	# Define functions for least squares fitting
	def ErfFitFunction(x,p):
		return (p[0]/2)*(1-sp.erf((x-p[1])/(np.sqrt(2)*p[2]))) #See Siegman, IEEE JOURNAL OF QUANTUM ELECTRONICS, VOL. 21, NO. 4, APRIL 1991
	
	def residuals(p,y,x):
		return y - ErfFitFunction(x,p)

	# Initialize vectors before the for loop
	waistsize = [0 for x in range(numfiles)] #[mm] Initialize waist size vector
	peakpower = [0 for x in range(numfiles)] #[W] Initialize laser beam power vector
	for i in range(numfiles):
		## Load Data for full knife edge measurement
		xdata,ydata = [], []
		filename = directoryname + filenames[i]
		inputfile = open(filename, 'r') # Open file for reading
		for line in inputfile:
			line = line.strip() # Read in a line at a time
			if not line.startswith("%"): # Only read a line if it doesn't start with "%", the comment symbol used in Matlab.
				columns = line.split()
				xdata.extend([float(columns[0])]) #[mm] Transverse position as knife edge cuts across beam
				ydata.extend([float(columns[1])]) #[mV] Power in beam in mV
		inputfile.close()

		ydata = np.abs(np.array(ydata))
		
		## Fitting
		# Guesses for data
		Aguess,cguess,oguess = (max(ydata)-min(ydata))/2, np.mean(xdata), np.mean(xdata)/(max(ydata)-min(ydata))
		#Non-linear least squares fit, passing x, y, and initial guesses
		plsq = optimize.leastsq(residuals,[Aguess,cguess,oguess],args = (ydata,xdata))
		bestFitValues = plsq[0]
		
		#Best fit values of amplitude (power), position of center of beam, and 2/e^2 radius of gaussian beam
		A,c,sigma = bestFitValues[0],bestFitValues[1],bestFitValues[2]
		
		# Vector of independent variable values to pass to the fit function
		xFit = np.arange(min(xdata)-1,max(xdata)+1,0.01)
		gaussianFromFit = (A/2)*np.exp(-(xFit-c)**2/(2*sigma**2)) # Vector of gaussian shape of cloud as function of position

		# Plotting
		plt.figure(i)
		plt.plot(xdata,ydata,'or') # Plot data
		plt.plot(xFit,ErfFitFunction(xFit,plsq[0]),'--b')
		plt.plot(xFit,gaussianFromFit,':g')

		plt.xlabel('Knife Edge Position [mm]')
		plt.ylabel('Measured Power [mV]')
		plt.text(np.mean(xdata)+0.2,np.mean(ydata),'1/e^2 radius = '+ str(abs(2*sigma)) + ' mm')
		plt.title('Beam Size at z = ' + str(measurementposition[i]) + ' cm')

		waistsize[i] = abs(2*sigma) #[mm] 1/e^2 radius of laser beam in mm
		peakpower[i] = abs(A) 		#[W] power of laser beam in watts

	
## Write data to output file
outputfilepath = directoryname + outputfilename; # Full path for output file
np.savetxt(outputfilepath, (measurementposition,waistsize),fmt='%1.4e') #[cm mm] Data is written to text file 