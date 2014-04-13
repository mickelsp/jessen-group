## Import Beam Profile Data From Beam Gage Software
# This script imports ASCII data exported from the Beam Gage software provided by Ophir Spiricon
# Data is in comma-delimited (.csv) format
# Then data is fit using a 2D gaussian function to extract x and y beam sizes

import matplotlib.pyplot as plt   # Import plotting functions from Matplotlib
import numpy as np                # Import numpy functions
import scipy.optimize as optimize # Import optimize function for least squares fitting

## Close and Clear Things
plt.close()
plotOn = 0    # 1 for plots to display; 0 to suppress plots

## Import Data
fileDirectory = “EnterDirectoryPath”Here
fileBatch = “EnterFileNameAndExtensionHere"

# Initially, just read list of file names. Line 30 actually imports data for each file
fileNames = []
zPosition = []
fo = open(fileDirectory+fileBatch)
values = fo.readlines()
for line in values:                 # First 16 lines are header; final line at 2065 is footer
  line = line.strip()               # Read in a line at a time
  if not line.startswith("%"): 	    # Only read a line if it doesn't start with "%", the comment symbol used in Matlab.
    columns = line.split()
    fileNames.append(columns[0])    # Append value of file name from this line to the list of file names
    zPosition.extend([float(columns[1])]) # Append value of z position to the list of z positions

fo.close()

## Loop Through all Files
numFiles = len(fileNames)
for i in range(0,numFiles):
  fileBase = fileNames[i]      # Select desired file amd then import its data
  beamData = np.loadtxt(fileDirectory+fileBase+'.csv',delimiter=',',usecols=range(1600)) # usecols command there b/c each row has a trailing comma that otherwise causes problems
  
  outputFileName = fileBase+'_FIT.txt'

  ## Transform Data
  beamData = np.flipud(beamData)    # Reorder the rows of the array because Python defined the origin differently than the Beam Gage software
  beamData = np.rot90(beamData,k=1) # Rotate the array by 90 degrees to make rows be the horizontal (x) data and columns the vertical (y) data
  beamData = beamData.ravel()

  ## Create X and Y Position Vectors
  pixelSizeCal = 4.4                                      #[microns/pixel] Size of each pixel
  numPixelsX,numPixelsY = 1200,1600                       #[pixels] Number of pixels in x and y directions
  positionX = np.arange(0,numPixelsX)*pixelSizeCal        #[microns] Position in horizontal direction
  positionY = np.arange(0,numPixelsY)*pixelSizeCal        #[microns] Position in vertical direction
  positionX,positionY = np.meshgrid(positionX,positionY)  #[microns] Create a grid with same dimensions as x and y axes
  
  indVars = (positionX,positionY)
  
  ## Fitting Data With a 2D Gaussian
  def Gaussian2D(indVars,amplitude,xPeak,yPeak,xSize,ySize,offset):
    x,y = indVars[0],indVars[1]
    zOut = amplitude * np.exp(-((x-xPeak)**2 / (2*xSize**2))-((y-yPeak)**2 / (2*ySize**2))) + offset
    return zOut.ravel()
    
  amplitudeGuess,xPeakGuess,yPeakGuess,xSizeGuess,ySizeGuess,offsetGuess = 1000,2600,3500,600,4000,10   # Initial guesses
  initialGuesses = (amplitudeGuess,xPeakGuess,yPeakGuess,xSizeGuess,ySizeGuess,offsetGuess)         # Initial guesses compiled
  
  popt, pcov = optimize.curve_fit(Gaussian2D,indVars,beamData,p0=initialGuesses)                    # Curve fitting
  
  beamFitted = Gaussian2D(indVars, *popt)                                                           # Create fitted data matrix
  
  ## Plot Beam Size Data
  if plotOn==1:
    plt.figure(i,figsize=(6,8))
    plt.imshow(beamData.reshape(numPixelsY, numPixelsX))    # Plot the raw data
    plt.contour(positionX[1,:]/pixelSizeCal,positionY[:,1]/pixelSizeCal,beamFitted.reshape(numPixelsY,numPixelsX),8,colors='w') # Super-impose fit values
    plt.xlabel('x [$\mu$m]')
    plt.ylabel('y [$\mu$m]')

  ## Write data to output file
  outputFilePath = fileDirectory + outputFileName     # Full path for output file
  np.savetxt(outputFilePath, (zPosition[i],popt[0],popt[1],popt[2],popt[3],popt[4],popt[5]),fmt='%1.4e')   #[arb um um um um arb] units