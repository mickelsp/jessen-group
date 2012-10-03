%%%Fitting function for beam size as a function of wavelength and inverse
%%%beam parameter
function w1=BeamSizeFunction(lambda,inverse)

%Calculation
w1=sqrt(-lambda./(pi.*imag(inverse))); %1/e^2 beam radius