%IMPLICITDISPERSE = a way of modelling dispersal implicitly
%   The method we used in Janssen et al., 2009 to model dispersion. 
%   For Rietkerks model this method is very efficient!
%   Because it uses fft use always a power of 2 for the number of cells (64x64 or 128x128 etc).
%
%Usage:
%   Include in model
%    X(1:64,1:64)' = f(X) + implicitdisperse(X,D)
%    or:
%    X(1:64,1:64)' = f(X) + implicitdisperse(X+b*Y,D)
%
%    See also model

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:27 $
function dA = implicitdisperse(ivar, avar, D, N)
%
%i_mmodel calls:
% implicitdisperse(ivar,"avar",D)
% and replaces "implicitdisperse(X,D)" in the model with:
% implicitdisperse(ivar,"X")
%
% Before each run i_run calls:
% implicitdisperse; = updating
%
%first argument should be the number of the state variable (for efficiency)
global g_grind;
if (nargin == 2)
   if isempty(g_grind.implicdisp{ivar}.D2s)
      dA=zeros(size(avar));
   else
      % Evolve in time t->t+dt
      a  = fftshift(fft2(avar));
      %implicit part of integration step
      aa = reshape(g_grind.implicdisp{ivar}.D2s * a(:), g_grind.implicdisp{ivar}.N, g_grind.implicdisp{ivar}.N);
      dA = (real(ifft2(ifftshift(aa))) - avar)./g_grind.solver.opt.MaxStep; 
      %divide by dt as dA is later multiplied with dt
   end;
elseif nargin == 3 %initiation
   %define dispersing state variable
   g_grind.implicdisp{ivar}.Name = avar;
   g_grind.implicdisp{ivar}.D = D;
   g_grind.implicdisp{ivar}.N = [];
   g_grind.implicdisp{ivar}.D2s = [];
elseif (nargin == 0) && isfield(g_grind,'implicdisp')
   % used internally to update g_grind, returns function in odefile
   % Wave numbers
   for ivar = 1:length(g_grind.implicdisp)
      if ~isempty(g_grind)&&~isempty(g_grind.implicdisp{ivar})
         dt=g_grind.solver.opt.MaxStep;
         if isempty(dt)
            error('GRIND:implicitdisperse:EulerNeeded','Use Euler integration for implicit dispersion');
         end;
         var = evalin('base',g_grind.implicdisp{ivar}.Name);
         N = length(var);
         log2N=log(N)/log(2);
         if round(log2N)~=log2N
            warning('GRIND:implicitdispers:power2','implicitdisperse is more efficient if the number of cells is a power of 2!');
         end;
         D = evalin('base',g_grind.implicdisp{ivar}.D);
         q = zeros(N, 1);
         for k = 1:N
            q(k) = (2 * (k - 1) / N - 1) * pi * N;
         end
         %Matrix for implicit part of integration, stored for later use
         Qmat = repmat(q.^2, 1, N) + repmat((q.^2)', N, 1);  %replicate and tile
         Qmat2 = spdiags(Qmat(:), 0, N^2, N^2);
         g_grind.implicdisp{ivar}.D2s = inv(speye(N^2) + dt * D * Qmat2);
         g_grind.implicdisp{ivar}.N = N;
      end;
   end;
end;
