function Z=myfitfunction(coeffs,data)
x = data; % get independent variable

%Initial Guesses
amplitude = coeffs(1);
center = coeffs(2);
width = coeffs(3);

%Fit Function
Z=amplitude.*exp(-2.*(x-center).^2./(width.^2));
