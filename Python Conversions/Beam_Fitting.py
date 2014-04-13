## Gaussian Beam Fitting
# For list of 1/e^2 sizes and positions, fit to beam propagation equation to get Rayleigh range, waist position, and waist size of beam
# Future additions: simultaneous fitting in both x and y transverse dimensions

import numpy as np                # Import numpy functions
import scipy.optimize as optimize # Import optimize function for least squares fitting
import matplotlib.pyplot as plt          # Import matplotlib functions

plt.close()

## Initialization
fontsize = 16         # Font size for plotting
laserlambda = 894e-9  #[m] wavelength of the laser beam
errorflag = 0         #1 to include error bars as weights in fitting; 0 otherwise

## Import data
# Data should be in 2-column format: longitudinal position (z) and 1/e^2 radius of the beam at that longitudinal position
directoryName = '/Users/work/Documents/Arizona/Nanofibers/Laser Profiling/Probe Beam Profiling/Before Fiber Link/'
inputFileName = 'ProbeLaser_FiberLink_ReverseOutput_20131008.txt'
fileName = directoryName + inputFileName
importedData = np.loadtxt(fileName)

zData = importedData[0]         #[cm] Longitudinal position along laser beam
beamSize = importedData[1]      #[mm] Transverse position across laser beam

## Data preparation
zData = zData*1e-2              #[m] converts z from cm to m
beamSize = beamSize*1e-3        #[m] converts beamSize from mm to m

## Fitting
#Fit to beam propagation equation
def GaussianPropagate(z,p):
  return p[0] * np.sqrt(1 + ( (z - p[1]) * laserlambda / (np.pi * p[0]**2) )**2)
  
def residuals(p,y,z):
  return y - GaussianPropagate(z,p)

waistGuess = 0.1e-3     #[m] Initial guess for waist value
zGuess = -0.7          #[m] Initial guess for position of waist

plsq = optimize.leastsq(residuals,[waistGuess,zGuess],args = (beamSize,zData))
fitWaist = plsq[0][0]         #[m] waist size; minimum beam size
fitWaistPosition = plsq[0][1] #[m] longitudinal position of minimum beam size

zFit = np.arange(-1,2,0.001)
beamPropagate = GaussianPropagate(zFit,[fitWaist,fitWaistPosition])#+0.00065
#beamPropagateGuess = GaussianPropagate(zFit,[0.35e-3,0.7])+0.00145

## Plot waist size versus distance from lens
plt.figure(1)
plt.plot(zData,beamSize,'sr',label="Measurements")
plt.plot(zFit,beamPropagate,'--b',label="Best Fit")
#plt.plot(zFit,beamPropagateGuess,'--g')
plt.xlabel('Longitudinal Position [m]')
plt.ylabel('Waist Size [mm]')
plt.text(0.05*np.mean(zData),1.1*np.mean(beamPropagate),'1/$e^2$ radius = '+ str(round(fitWaist,5)) + ' m')
plt.text(0.05*np.mean(zData),1.05*np.mean(beamPropagate),'z$_0$ = '+ str(round(fitWaistPosition,4)) + ' m')
plt.legend(bbox_to_anchor=(0, 0, 1, 1), bbox_transform=plt.gcf().transFigure)
plt.title('Beam Profile of Nd:YAG Laser')
