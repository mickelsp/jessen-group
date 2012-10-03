function Z=BosonFugacityFitFunction(coeffs,data)
peakOD=coeffs(1);
sigmax=coeffs(2);
sigmay=coeffs(3);
xoffset=coeffs(4);
yoffset=coeffs(5);
offset = coeffs(6); 
slopex=coeffs(7);
slopey=coeffs(8);
%fugacity=coeffs(9);
fugacity=abs(coeffs(9));
%fugacity=1;

accuracy= .01;
x=data(:,1);  %Data split into x and y vectors
y=data(:,2);
w=data(:,3);  %Extra vector for weighting if needed
numberofpoints=data(:,4);%matrix in which each entry is the number of points 
sqrtnumbpoints=sqrt(numberofpoints(1,1)); %sqrt of number of points

%Z = (offset+slopex*(x-xoffset)+slopey*(y-yoffset)+(peakOD*polylog(2,abs(sin(fugacity*(pi/2))).*exp(-(x-xoffset).^2/(2*sigmax^2)).*exp(-(y-yoffset).^2/(2*sigmay^2)),accuracy)/polylog(2,abs(sin(fugacity*(pi/2))),accuracy)))./(w*sqrtnumbpoints)
Z = (offset+slopex*(x-xoffset)+slopey*(y-yoffset)+(peakOD*polylog(2,fugacity.*exp(-(x-xoffset).^2/(2*sigmax^2)).*exp(-(y-yoffset).^2/(2*sigmay^2)),accuracy)/polylog(2,fugacity,accuracy)))./(w*sqrtnumbpoints);
%Z = (offset+slopex*(x-xoffset)+slopey*(y-yoffset)+(peakOD*dilog(fugacity.*exp(-(x-xoffset).^2/(2*sigmax^2)).*exp(-(y-yoffset).^2/(2*sigmay^2)))/dilog(fugacity)))./(w*sqrtnumbpoints)

end

