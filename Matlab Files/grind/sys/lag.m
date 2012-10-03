%LAG   Time lag in a differential equation.
%   Use a time lag for some state variable (e.g. May 1976). Use
%   this function in differential equations. The delay differential
%   equation (DDE) is solved with the dde23 solver.
%   You can use the same function among others in time plots (see out) to plot 
%   lagged state variables, auxiliary variables or functions.
% 
%   Usage:
%   LAG('statevar',TIMELAG)
%   statevar = name of state variable.
%   TIMELAG   = lagging time (time units).
%
%   See also model, dde23, out

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:27 $

% this function is replaced in the ODE function as it is produced by DDE
% it can however be used afterwards as function, var can then also be a function
function res=lag(var,timelag)
global g_t;
if ischar(var)
   var=i_getoutfun(var);
end;
if size(var,1)~=size(g_t)
   error('GRIND:lag:ArgError','Error in lag function, var should be a string or a matrix of size g_t');
end;
res=interp1(g_t,var,g_t-timelag);


