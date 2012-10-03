%FFTCONV   Convolution using fft (for Integro-Difference Equations)
%   This function implements convolution of a dispersion kernell using a fft algorithm.
%   Use this function to implement dispersion in an integro-difference equation (Powell, xxxx).
%   It has to be used in combination with setevent in a 1D or 2D spatial explicit model.
%
%   Usage:
%   first use: fftconv(X,L,pars,distribution) in which X is the matrix, L is the size of
%   each grid cell (can be two dimensional [LX,LY]), pars is a vector with parameters for 
%   each distribution, distribution can be 'normal','double' or 'ring'.
%   other use: Xnew=fftconv(X);
%   dispersion kernells (see Etienne, et al. 2002):
%   'normal': random walk parameters: pars=[D*deltat]
%   'double': double Normal distribution: pars = [D*deltat]
%   'ring': ring random distribution: pars = [D*deltat,rho]
%   
%
%   Examples:
%   Add this to a model (parameter panel) to immplement an Integro difference equation.
%   setevent('simpleevent',0,'fftconv(X,L,[D*deltat],''normal'');',NaN);
%   setevent('simpleevent',deltat,'X=fftconv(X);',deltat);
%
%
%  See also setevent, model

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:26 $
function [P1] = fftconv(P, L, pars,kernell)
global g_grind;
if ((nargin==4))||~isfield(g_grind,'fK')
   g_grind.fK = getdistrib(P, L, pars,kernell);
end;
if size(P, 2) <= 1
   P1 = real(fftshift(ifft(g_grind.fK .* fft(P))));
else
   P1 = real(fftshift(ifft2(g_grind.fK .* fft2(P))));
end;
   
function fK = getdistrib(P, L, pars,kernell)
[nx, ny] = size(P);
if ny <= 1
   if length(L) > 1
      L = L(1);
   end;
   xlr = L * nx / 2;
   x = linspace(-xlr, xlr-L, nx)';
   % add here 1D dispersion kernells
   if strcmpi(kernell,'normal') %pars(1) = D; pars(2) = deltat (or pars(1)=D*deltat)
      K = 1 ./ sqrt(4*pi*prod(pars)) .* exp(-x.^2  ./ (4 .* prod(pars)));
   elseif strcmpi(kernell,'double') %pars(1) = D; pars(2) = deltat (or pars(1)=D*deltat)
      sigma=sqrt(2*prod(pars)); 
      lambda=2*sqrt(2)/(sigma*sqrt(pi));
      K = lambda^2/(2*pi)* exp(-lambda*(x ));
   elseif strcmpi(kernell,'ring') %pars(1)=D*deltat pars(2)=rho
      K = exp(-(x-pars(2)).^2 ./ (4 .* pars(1)));
   else
      error('GRIND:fftconv:UnknownKernell','Dispersion kernell unknown');
   end;
   K = K / (L * trapz(K)); %normalize;
   fK = L * fft(K);
else
   if length(L) <= 1
      L(2) = L(1);
   end;
   xlr = L(1) * nx / 2;
   x1 = linspace(-xlr, xlr-L(1), nx)';
   ylr = L(2) * ny / 2;
   y1 = linspace(-ylr, ylr-L(2), ny)';
   [x, y] = meshgrid(x1, y1);
   % add here 2D dispersion kernells (see Etienne 2002)
   if strcmpi(kernell,'normal') %pars(1) = D; pars(2) = deltat (or pars(1)=D*deltat)
      K = 1 ./ (4*pi*prod(pars)) .* exp(-(x.^2 +  y.^2) ./ (4 .* prod(pars)));
   elseif strcmpi(kernell,'double') %pars(1) = D; pars(2) = deltat (or pars(1)=D*deltat)
      sigma=sqrt(2*prod(pars)); 
      lambda=2*sqrt(2)/(sigma*sqrt(pi));
      K = lambda^2/(2*pi)* exp(-lambda*(sqrt(x.^2 +  y.^2)));
   elseif strcmpi(kernell,'ring') %pars(1)=D*deltat pars(2)=rho
      K = exp(-(sqrt(x.^2 +  y.^2)-pars(2)).^2 ./ (4 .* pars(1)));
   else
      error('GRIND:fftconv:UnknownKernell','Dispersion kernell unknown');
   end;
   K = K / (prod(L) * trapz(trapz(K))); %normalize
   fK = prod(L) * fft2(K);
end


