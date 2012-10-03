%DWIENER   Add stochasticity to a differential equation (Wiener process)
%  Consider the stochastic differential equation:
%  dX = f(X) dt + g(X) dW
%  in which 
%     X is the state variable 
%     f(X) the "drift term" as function of X
%     g(X) the "diffusive term" as function of X (or a constant in case of additive noise)
%     dW the increment of the Wiener process or Brownian motion
%  
%  This model can be defined in GRIND as:
%  X' = f(X) + dwiener(g(X),g'(X))
%  in which g'(X) is the derivative of g(X) to X. (thus zero for additive noise)
%  If this function is added to a differential equation, the solver is set to Euler integration
%  and the model is solved with the explicit Milstein method (if g' is defined) or with the Euler-Maruyama
%  scheme for additive noise. (both are based on the Ito integral).
%
%
%  Usage:
%  DWIENER(SIGMA) - adding additive noise with Normal distribution with standard deviation SIGMA.
%  DWIENER(GX,GXACCENT) - adding multiplicative noise solved with Milstein method.
%  DWIENER(GX) - if GXACCENT is ignored the less efficient Euler-Maruyama method is used.
%
%  See also modelpanel, model, rednoise
%  

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:26 $
function res=dwiener(gY,dgY_dY,at)
global g_grind t;
if nargin==2
   at=dgY_dY;
   dgY_dY=0;
end;
if ~strcmpi(g_grind.solver.name,'euler')
   if g_grind.solver.isdiffer
      error('GRIND:dwiener:diffEquation','dwiener cannot be used for a difference equation, use randn() instead');
   end;
   error('GRIND:dwiener:NoEuler','dwiener process needs Euler integration');
end;
if length(at)>1
   %remake a similar (BUT NOT THE SAME!) data set if used in a function.
   gY=gY';
   dgY_dY=dgY_dY';
   h=1;
   dW=randn(length(at),length(gY));
   warning('GRIND:dwiener:newseries','Note that dwiener generated a new (but similar) noise series');
else
   h=g_grind.solver.opt.MaxStep;
   if at== t + 0.00001 %null uses this
      dW=zeros(size(gY));
   else
      dW=randn(size(gY)).*sqrt(h);
   end;
end;
%Milstein scheme (reduces to Euler-Maruyama if d(gY)/dY=0)
res=(gY.*dW+gY.*dgY_dY.*(dW.^2-h))./h;% divide by h as the Euler routine multiplies with h
