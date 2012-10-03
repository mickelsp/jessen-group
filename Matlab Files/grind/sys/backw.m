%BACKW   Simulate the model backwards
%   Run the model with negated right-hand sides (=backwards). 
%   For difference equations this is more complicated due to the 
%   discrete steps. Therefore a local optimizer is used to find 
%   the next step back in time.
%   This way you can find unstable equilibria (see also: findeq) or a separatrix.
%   (see: perturb).
%
%   Usage:
%   BACKW - run backwards with g_grind.ndays time steps (see also SIMTIME for setting the simulation duration).
%   BACKW N - run backwards for N steps
%
%   See also ru, simtime, solver, perturb, addmode

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:26 $
function backw(ndays,dotrunc)
global g_grind;
i_parcheck;
if nargin == 0
   ndays = g_grind.ndays;
else
   ndays = i_checkstr(ndays);
end;
if nargin<2
   dotrunc=1;
end;
g_grind.solver.backwards =  ~g_grind.solver.backwards;
oldtrunc=g_grind.truncate;
g_grind.truncate=dotrunc;
try
   ru(ndays)
   g_grind.solver.backwards = ~g_grind.solver.backwards;
   g_grind.truncate=oldtrunc;
catch err
%   err=lasterror;
   g_grind.solver.backwards = ~g_grind.solver.backwards;
   g_grind.truncate=oldtrunc;
   rethrow(err);
end;
