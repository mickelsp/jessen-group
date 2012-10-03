%%%Fitting function for knife edge measurements (fit to error function)
function Z=myfitfunction(coeffs,data)
x = data(:, 1); % get independent variable; power measurements

%Initial Guesses
amplitude1=coeffs(1); %total power in the beam
center1 = coeffs(2); %position of center of beam
w1=coeffs(3); %1/e^2 beam radius

%Fit Function
Z=(amplitude1./2).*(1-erf((x-center1)./(sqrt(2).*w1))); %See Siegman, IEEE JOURNAL OF QUANTUM ELECTRONICS, VOL. 21, NO. 4, APRIL 1991