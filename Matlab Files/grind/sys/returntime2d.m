%RETURNTIME2D   Plot of the time necessary to reach a stable node
%   Create a 2D contourplot of the phase plane containing the time to reach any stabe equilibrium. 
%   
%   Usage:
%   RETURNTIME2d - estimates the returntime based on the current SIMTIME settings and a default
%   value for the maximum change in equilibrium (1E-8). 
%   [V,X,Y]=RETURNTIME2D write the results also to V X and Y for further use for instance in SURF or
%   other MATLAB functions.
%   RETURNTIME2D ERR - use the maximum change of ERR.
%   RETURNTIME2D ERR MAXT - use ERR and simulate MAXT timesteps to find an equilibrium.
%   
%   
%   See also returntime, ru, simtime, findeq

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:28 $
function [Vect,X,Y]=returntime2d(err,npoints)
global g_grind;
i_parcheck;
if nargin ==0
   err = 1E-8;
else
   err =  i_checkstr(err);
end;
if nargin < 2
   npoints = 50;
else
   npoints = i_checkstr(npoints);
end;
iX = i_getno(g_grind.xaxis.var);
iY = i_getno(g_grind.yaxis.var);
if (isempty(iX.no) || isempty(iY.no))
   ax ?;
   errordlg('Cannot returntime plot if there are no state variables/parameters on the axes.');
   error('GRIND:returntime2d:InvalidAxes','Cannot returntime plot if there are no state variables/parameters on the axes');
end
oldX = evalin('base',g_grind.xaxis.var);
oldY = evalin('base',g_grind.yaxis.var);
try
   Vect = zeros(npoints, npoints);
   X = zeros(npoints, npoints);
   Y = zeros(npoints, npoints);
%   t = 0;
   minX=g_grind.xaxis.lim(1);
   if minX<0.0001, minX=0.0001; end;
   maxX=g_grind.xaxis.lim(2);
   minY=g_grind.yaxis.lim(1);
   if minY<0.0001, minY=0.0001; end;
   maxY=g_grind.yaxis.lim(2);
   incrY=(maxY - minY) / (npoints-1);
   incrX=(maxX - minX) / (npoints-1);
%   isdiffer = g_grind.solver.isdiffer;
   for y1 = 1:npoints
      py = (y1 - 1) * incrY + minY;
      assignin('base',g_grind.yaxis.var,py)
      for x1 = 1:npoints
         px = (x1 - 1) * incrX + minX;
         assignin('base', g_grind.xaxis.var, px)
         Vect(x1, y1) = returntime(err);
         X(x1, y1) = px;
         Y(x1, y1) = py;
      end;
   end;
   assignin('base', g_grind.xaxis.var, oldX)
   assignin('base', g_grind.yaxis.var, oldY)
catch err
%   err=lasterror;
   assignin('base', g_grind.xaxis.var, oldX)
   assignin('base', g_grind.yaxis.var, oldY)
   rethrow(err);
end;
[H,new] = i_makefig('phase2');
if new
set(H, 'WindowButtonDown', 'i_callb(''mdown'')');
set(H, 'WindowButtonMotionFcn', 'i_callb(''mmove'')');
end;
set(H, 'Name', 'Phase plane');
oldhold = ishold;
hold on;
surf(X, Y, Vect);
shading flat; 
set(gca,'XLim', g_grind.xaxis.lim);
set(gca,'YLim', g_grind.yaxis.lim);
if ~oldhold
   hold off;
end;



