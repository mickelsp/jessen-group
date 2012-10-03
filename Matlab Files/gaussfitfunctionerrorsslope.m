function Z = myfitfunction(coeffs, data);
amplitude = coeffs(1);
sigx = coeffs(2);
sigy = coeffs(3);
xoffset=coeffs(4);
yoffset=coeffs(5);
offset = coeffs(6); % Extract the coefficients from the vector
slopex=coeffs(7);
slopey=coeffs(8);

x = data(:, 1); % Split the data matrix into x and y vectors
y = data(:, 2);
w = data(:, 3); %matrix of errors
numberofpoints=data(:,4); %matrix in which each entry is the number of points 
sqrtnumbpoints=sqrt(numberofpoints(1,1)); %sqrt of number of points
%offset=1000;
%amplitude=-5000;
%xoffset=42;
%yoffset=46;
%slopex=-.0001;
%slopey=.0003;

Z =((offset+slopex*(x-xoffset)+slopey*(y-yoffset)+amplitude*exp(-(((x-xoffset).^2)/(2*sigx^2))-(((y-yoffset).^2)/(2*sigy^2))))./(w*sqrtnumbpoints)); %We're trying to fit z = f(x, y) so compute f(x, y) 