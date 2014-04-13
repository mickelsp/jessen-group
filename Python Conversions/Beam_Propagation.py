## Gaussian Beam Propagation
# For some beam size and wavelength of light, determine the beam's propagation characteristics. If there is a lens, predict the changed path.
# Future additions: propagation of astigmatic beams (x and y waist locations different)

import numpy as np # Import numpy functions
import scipy.optimize as optimize # Import optimize function for least squares fitting
import matplotlib.pyplot as plt # Import matplotlib functions

plt.close()

## Initialization
laserlambda = 1064e-9 #[m] Wavelength of the beam
z = np.arange(-1,2,0.001) #[m] Longitudinal position

## ABCD Matrix Functions for Propagation
def ABCDmatrix_FreeSpace(inverse1,z):
    inverse2 = inverse1/(1+z*inverse1) # New inverse beam parameter after distance, z
    return inverse2
  
def ABCDmatrix_ThinLens(inverse1,f):
    inverse2 = -1/f + inverse1 # New inverse beam parameter after thin lens of focal length, f
    return inverse2
    
def BeamSizeFunction(wavelength,inverse1):
    w1 = np.sqrt(-laserlambda/(np.pi*inverse1.imag)) # 1/e^2 beam radius as a function of wavelength and the inverse beam parameter
    return w1

## Incoming Beam Characteristics
# Incoming beam characteristics could be extracted from measurement or they could be arbitrarily chosen if one were designing a system from scratch.
w0 = 3.5e-4 #[m] 1/e^2 radius of beam at its maximum intensity
z0 = 0.6899 #[m] position of the minimum beam size in lab coordinate system
zR = np.pi*w0**2/laserlambda #[m] Rayleigh range of the beam

qinitialinverse = -1J*laserlambda/(np.pi*w0**2) #[1/m] Value of the inverse of the complex beam parameter at the waist of the incoming beam
qinitialinversevector = ABCDmatrix_FreeSpace(qinitialinverse,z-z0) #[1/m] Vector of inverse of the complex beam parameter as a function of position for incoming beam

incomingbeamsize = BeamSizeFunction(laserlambda,qinitialinversevector) #[m] Beam size as a function of z position

## Lens 1 Beam Characteristics
f1 = 150e-3 #[m] Focal length of lens 1
z1 = 0.010 #[m] Position of lens 1 relative to waist of incoming beam
zlens1 = z0+z1 #[m] Position of lens 1 relative to lab coordinate system

qbeforelens1inverse=ABCDmatrix_FreeSpace(qinitialinverse,z1) #[1/m] Value of the inverse of the complex beam parameter just before lens 1
qafterlens1inverse = ABCDmatrix_ThinLens(qbeforelens1inverse,f1) #[1/m] Value of the inverse of the complex beam parameter just after lens 1
qafterlens1inversevector = ABCDmatrix_FreeSpace(qafterlens1inverse,z-zlens1) #[1/m] Vector of the inverse of the complex beam parameter after lens 1

afterlens1beamsize = BeamSizeFunction(laserlambda,qafterlens1inversevector) #[m] 1/e^2 beam radius as a function of z position

## Lens 2 Beam Characteristics
f2 = 38.1e-3 #[m] Focal length of lens 2
z2 = 0.187 #[m] Position of lens 2 relative to lens 1 waist
zlens2 = zlens1 + z2 #[m] Position of lens 2 relative to lab coordinate system

qbeforelens2inverse=ABCDmatrix_FreeSpace(qafterlens1inverse,z2) #[1/m] Value of the inverse of the complex beam parameter just before lens 2
qafterlens2inverse = ABCDmatrix_ThinLens(qbeforelens2inverse,f2) #[1/m] Value of the inverse of the complex beam parameter just after lens 2
qafterlens2inversevector = ABCDmatrix_FreeSpace(qafterlens2inverse,z-zlens2) #[1/m] Vector of the inverse of the complex beam parameter after lens 2

afterlens2beamsize = BeamSizeFunction(laserlambda,qafterlens2inversevector) #[m] 1/e^2 beam radius as a function of z position

## Plot Beam Propagation
zeropoint=300 # Position at z vector where zero lives
plt.figure(1)
plt.plot(z*1e2,incomingbeamsize*1e3,'-k',linewidth=1.5,label="Input Beam") # Plots size of beam after lens 2 over z range
plt.plot(z*1e2,afterlens1beamsize*1e3,'-g',linewidth=1.5,label="After Lens 1") # Plots size of beam after lens 3 over z range
plt.plot(z*1e2,afterlens2beamsize*1e3,'-r',linewidth=1.5,label="After Lens 2") # Plots size of beam after lens 3 over z range
plt.xlabel('Longitudinal Position [cm]')
plt.ylabel('1/$e^2$ Beam Radius [mm]')
# plt.xlim((0,40))
# plt.ylim((0,4))
plt.legend(bbox_to_anchor=(0, 0, 1, 1), bbox_transform=plt.gcf().transFigure)
plt.title("Nd:YAG Output Beam Profile")