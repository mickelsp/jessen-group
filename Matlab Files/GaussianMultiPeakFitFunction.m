function Z=myfitfunction(coeffs,data)
x=data(:, 1); % get independent variable
y=data(:, 2); % 
width=data(1,4); %get width constant

%Initial Guesses
center1=coeffs(1);
center2=coeffs(2);
center3=coeffs(3);
%width=coeffs(4);
amplitude=coeffs(4);

%width=2;
%amplitude=25;

%Fit Function
Z=amplitude*1./(1+(2*(x-center1)/width).^2)+amplitude*1./(1+(2*(x-center2)/width).^2)+amplitude*1./(1+(2*(x-center3)/width).^2);
