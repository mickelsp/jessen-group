global g_Y;
n=65536;
Y = fft(g_Y,n);

%The power spectral density, a measurement of the energy at various frequencies, is

Pyy = Y.* conj(Y) / n;

%Graph the first 257 points (the other 255 points are redundant) on a meaningful frequency axis. 

f = 1000*(0:n/2)/n;
hold on;
plot(f,Pyy(1:n/2+1,:));
set(gca,'yscale','log','xscale','log');



