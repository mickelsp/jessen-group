### Gaussian Beam Fitting
#For list of 1/e^2 sizes and positions, fit to beam propagation equation
#to get Rayleigh range, waist position, and waist size of beam
#Future additions: simultaneous fitting in both x and y transverse dimensions

import numpy as np #import numpy functions
import scipy as sp #import scipy functions
#import matplotlib as mpl #import matplotlib functions
import csv #import csv functions

###Initialization
fontsize = 16;
laserlambda = 690e-9; #[m] wavelength of the laser beam
errorflag = 1; #1 to include error bars as weights in fitting; 0 otherwise

###Import data
#Data should be in 2 column format: longitudinal position (z) and 1/e^2 radius of the beam at that longitudinal position
directoryname = '/Users/work/Documents/Nanofibers/Laser Profiling/Probe Beam Profiling/After Nanofiber/'; 
inputfilename = 'ProbeLaser_NanofiberOutput_Telescope1_20130906.dat';
filename = directoryname + inputfilename;
z = []; beamsize=[]; #initialize vectors
inputfile = open(filename, 'r') #open file for reading
for line in inputfile:
	line = line.strip() #read in a line at a time
	if not line.startswith("%"): #only read a line if it doesn't start with "%", the comment symbol I've been using in Matlab.
		columns = line.split()
		z = z + [float(columns[0])]; #append value of z from this line to the list of z values
		beamsize = beamsize + [float(columns[1])];
inputfile.close()

###Convert lists into numeric vectors
z = np.array(z); #[cm]
beamsize = np.array(beamsize); #[um]

###Data preparation
z = z*1e-2; #[m] converts z from cm to m
beamsize = beamsize*1e-3; #[m] converts beamsize from mm to m
lambdavector = np.zeros(len(z));
lambdavector[0] = laserlambda;

###Fitting
#Fit to beam propagation equation
def func(z,waist,z0):
  return waist * sqrt(1 + ( (z - z0) * laserlambda / (pi * waist**2) )**2)
  
z = np.linspace(-4,0.001,4)
fitfunc = func(z,waist,zo);
fitfuncn = fitfunc + 0.2*np.random.normal(size=len(z));

popt, pcov = optimize.curve_fit(func,z,fitfuncn)

InitialGuess=np.array([1e-4,3]); #[m m] initial guess for waist value and position of waist
#datamatrix = [z(:) lambdavector(:)]; %[m m] data matrix contains both the values to be fitted and the fixed parameter (laser wavelength)
#if errorflag==1
#    errorvector = beamsize./beamsize; %weight all points equally in absence of better information
#    %errorvector=[0.8 0.4 0.4 0.2 0.1];
#    [P,r,J]=nlinfitweight(datamatrix,beamsize(:),@fittogaussianbeampropagation,InitialGuess,errorvector(:));
#else
#    [P,r,J]=nlinfit(datamatrix,beamsize(:),@fittogaussianbeampropagation,InitialGuess);
#end

#fitwaist = P(1) %[m] waist size; minimum beam size
#fitwaistposition = P(2) %[m] longitudinal position of minimum beam size

###Plot waist size versus distance from lens
#figure(1)
#plot(measurementposition(:),waistsize(:),'sr','MarkerFaceColor','r');
#hold on
#set(gca,'FontSize',fontsize,'FontWeight','bold');
#xlabel('Longitudinal Position [cm]','FontSize',fontsize,'FontWeight','bold');
#ylabel('Waist Size [mm]','FontSize',fontsize,'FontWeight','bold');
#%xlim([0 1.13]);
#ylim([0 max(waistsize(:))+0.1.*max(waistsize(:))]);
#%legend('Vertical Size','Horizontal Size','Location','best');
#%title('Beam Profile of 2988 Laser','FontSize',fontsize,'FontWeight','bold');

#Manually test
print z
print beamsize
print InitialGuess
print 'status = done'
