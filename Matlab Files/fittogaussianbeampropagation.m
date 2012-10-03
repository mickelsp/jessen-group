%% Gaussian beam propagation fitting function
function beamsize = fittogaussianbeampropagation(coeffs, data)
waist = coeffs(1); %waist is a fit parameter
z0 = coeffs(2); %waist position is a fit parameter

z = data(:,1); % Split the data matrix into x and y vectors
lambda = data(1,2); %[m] wavelength of light: constant, not a fitting parameter

beamsize = waist.*sqrt(1+((z-z0).*lambda./(pi.*waist.^2)).^2);