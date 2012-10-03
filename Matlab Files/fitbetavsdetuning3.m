function y = fitbetavsdetuning(coeff,inputmatrix)
amp=coeff(1);%fit is for scaling factor; best value of loptical from fitting the losses is amp.*lopticalfitting

delta=inputmatrix(:,1);
betavalues=inputmatrix(:,2);
intensityvalues=inputmatrix(:,3);
lopticalfitting=inputmatrix(1,4); %best value of loptical determined by fitting cross section data

h=6.626*10^(-34); %Planck's constant
hbar=h./(2.*pi); %divide by 2pi
massSr=88*1.672*10^(-27); %strontium mass
bohrradius=5.29e-11; %Bohr radius in meters
atomiclinewidth=7500; %7.5 kHz linewidth of intercombination line
gammanatural=2.*atomiclinewidth; %linewidth of PAS line; twice the atomic linewidth
deltae=-0.435e6; %binding energy of second shallowest PAS line in Hz (24 MHz).

%%%Fit function doesn't actually calculate beta.  We're fitting to minimize
%%%the difference between the function and the data (target value for this
%%%difference is zero).  This allows us to weight each point by its own
%%%value.
y = (4.*pi.*hbar./(massSr./2).*amp.*lopticalfitting.*bohrradius.*(gammanatural./(delta.*1e9-deltae)).^2-(betavalues./intensityvalues)).^2./(betavalues./intensityvalues).^2;