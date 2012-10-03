%STABIL   Stabilize
%   Run the current model without showing results and keep the final state as
%   initial value
% 
%   Usage:
%   STABIL - runs for 1000 days.
%   STABIL N - runs for N days.
% 
%   See also ru, ke

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:29 $
function stabil(ndays,silent)
global t g_grind;
if nargin == 0
   ndays = 1000;
else
   ndays = i_checkstr(ndays);
end;
if nargin<2
   silent=0;
end;
i_parcheck;
N0 = i_initvar;
oldstep = g_grind.tstep;
try
   if ~(g_grind.solver.isdiffer||g_grind.solver.haslag)
      g_grind.tstep = 2;
   else
      g_grind.tstep=NaN;
   end;
   g_grind.solver.opt.OutputFcn = [];
   i_ru(g_grind.odefile, t, ndays, N0, 1);
   g_grind.tstep = oldstep;
   ke;
   if ~silent
      fprintf('Simulated %d days.\n', ndays);
   end;   
catch err
%   err=lasterror;
   g_grind.tstep = oldstep;
   rethrow(err);
end;
