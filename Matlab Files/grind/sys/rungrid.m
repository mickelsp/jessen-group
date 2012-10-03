%RUNGRID   Create grid of trajectories
%   Generates a grid of trajectories in a 2D phase plane. These trajectories start at 
%   regular intervals in the state space.
%
%   Usage:
%   RUNGRID - creates a 5x5 grid and runs with the default number of days
%   RUNGRID N - creates a NxN grid
%   RUNGRID NX NY - creates a NX x NY grid. If one of these arguments are smaller than
%   1, it is not varied. Instead the current initial conditions are used.
%   RUNGRID NX NY NDAYS - runs for NDAYS.
%   RUNGRID NX NY NDAYS 1 - runs for NDAYS in two directions.
%
%   See also ru, phas, null 

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:28 $
function rungrid(nx, ny, days, twodir)
global g_grind;
i_parcheck;
if nargin == 0
   nx = 5;
end;
if nargin <= 1
   ny = nx;
end;
if nargin <= 2
   days = g_grind.ndays;
end;
if nargin<4
   twodir=0;
end;
nx = i_checkstr(nx);
ny = i_checkstr(ny);
days = i_checkstr(days);
oldndays = g_grind.ndays;
g_grind.ndays = days;
olddraw = g_grind.drawnow;
g_grind.drawnow = 0;
OldN0 = i_initvar;
N0 = OldN0;
ix = i_getno(g_grind.xaxis.var);
iy = i_getno(g_grind.yaxis.var);
if isempty(iy.no) && isempty(ix.no)
   errordlg('Can not set initial variables, because there are no state variables on the axes');
else
   if (nx > 1) && (ny > 1)
      for X = g_grind.xaxis.lim(1) + 0.01:(g_grind.xaxis.lim(2) - g_grind.xaxis.lim(1)) / (nx - 1):g_grind.xaxis.lim(2)
         for Y = g_grind.yaxis.lim(1) + 0.01:(g_grind.yaxis.lim(2) - g_grind.yaxis.lim(1)) / (ny - 1):g_grind.yaxis.lim(2)
            if ~DoRun(N0, ix, X, iy, Y, twodir);
               return;
            end;      
        end;
      end;
   elseif (nx <= 1) &&(ny>1)
      X = NaN;
      for Y = g_grind.yaxis.lim(1) + 0.01:(g_grind.yaxis.lim(2) - g_grind.yaxis.lim(1)) / (ny - 1):g_grind.yaxis.lim(2)
        if ~DoRun(N0, ix, X, iy, Y, twodir);
           return;
        end;      
      end;
   elseif (ny <= 1) &&(nx>1)
      Y = NaN;
      for X = g_grind.xaxis.lim(1) + 0.01:(g_grind.xaxis.lim(2) - g_grind.xaxis.lim(1)) / (nx - 1):g_grind.xaxis.lim(2)
         DoRun(N0, ix, X, iy, Y, twodir);
      end;
   else
      Y=NaN;
      X=NaN;
      DoRun(N0, ix, X, iy, Y, twodir);
   end;
end;
i_keep(OldN0);
g_grind.drawnow =  olddraw;
g_grind.ndays = oldndays;

function Ok=DoRun(N0, ix, X, iy, Y, twodir)
N0 = setaxisval(N0, ix, X);
N0 = setaxisval(N0, iy, Y);
i_keep(N0);
Ok=1;
ru;
if twodir
   backw;
end;
ud=get(gcf,'userdata');
if ~isempty(ud)&&isfield(ud,'stop')
   Ok=~ud.stop;
end;
drawnow;

function N0 = setaxisval(N0, ix, aval)
global g_grind;
if ~isempty(ix.no) && ~isnan(aval)
   if ix.ispar
      evalin('base', g_grind.pars{ix.no});
      assignin('base', g_grind.pars{ix.no}, aval);
   else
      N0(ix.no) = aval;
   end;
end;

