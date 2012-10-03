function Z=myfitfunction(coeffs,data)
x = data(:, 1); % get independent variable


%%%FOR THREE PEAKS:
amplitude1=coeffs(1);
center1 = coeffs(2);
width1=coeffs(3);

amplitude2=coeffs(4);
center2=coeffs(5);
width2=coeffs(6);

offset=coeffs(7);

Z=amplitude1*1./(1+(2*(x-center1)/width1).^2)+amplitude2*1./(1+(2*(x-center2)/width2).^2)+offset;


%%%FOR ONLY ONE PEAK:
% amplitude1=coeffs(1);
% center1 = coeffs(2);
% width1=coeffs(3);
% offset=coeffs(4);
% Z=amplitude1*1./(1+(2*(x-center1)/width1).^2)+offset;