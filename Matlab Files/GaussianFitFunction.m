function Z=myfitfunction(coeffs,data)
x = data(:, 1); % get independent variable

%Initial Guesses
amplitude1=coeffs(1);
center1 = coeffs(2);
width1=coeffs(3);
offset=coeffs(4);

%Fit Function
Z=amplitude1*1./(1+(2*(x-center1)/width1).^2)+offset;
