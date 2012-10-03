%PLOTDIFF   Plot 1D differential equation
%   Create a plot of a differential equation of a one dimensional
%   system. On the y-axis the growth is plotted.
%
%   Usage:
%   PLOTDIFF - plots the variable of the x axis with the range of the x-axis.
%   PLOTDIFF VAR1 LIM - plots the state variable VAR1 with a range of LIM. 
%
%   See also null, phas, ax, plotreldiff

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:28 $
function [aX,aY] = plotdiff(avar, alim)
global g_grind;
npoints=500;
i_parcheck;
if nargin == 0
   avar = g_grind.xaxis.var;
end;
if nargin <= 1
   alim = g_grind.xaxis.lim;
else
   alim = i_checkstr(alim);
end;
%if size(g_grind.statevars,2)>1
%   errordlg('Warning: This command is designed for 1D differential equations');
%end;
iX = i_varno(avar);
if isempty(iX)
   error('GRIND:plotdiff:NoStatevar','Can only create plot if state variable is on the axis');
end;
N0 = i_initvar;
X = (alim(1):(alim(2) - alim(1)) /npoints:alim(2));
Y = zeros(1, size(X, 2));
for i = 1:size(X, 2)
   N0(iX) = X(i);
   yy = feval(g_grind.odefile, 1, N0);
   Y(i) = yy(iX);
end;
if g_grind.solver.backwards
   Y = -Y;
end;
if nargout > 0
   aX = X;
   aY = Y;
else
   [H,new] = i_makefig('phase1');
   set(H, 'Name','Plot of 1D differential equation')
   if new
   set(H, 'WindowButtonDown', 'i_callb(''mdown'')');
   set(H, 'WindowButtonMotionFcn', 'i_callb(''mmove'')');
   end;
   oldhold = ishold;
   hold on;
   plot(X, Y);
   i_plotdefaults(H)
   xlabel(i_disptext(avar));
   ylabel([i_disptext(avar) '''']);
   plot(X, zeros(1, size(X, 2)), 'k');
   if g_grind.statevars.dim > 1
      title(i_disptext(['Valid for ' i_othervars(N0, iX)]));
   end;
   set(gca, 'XLim', alim);
   if ~oldhold
      hold off;
   end;
end;
if g_grind.solver.isdiffer
   i_warningdlg('GRIND:plotdiff:differenceeq','This function is designed for differential equations, please use <a href="itermap">itermap</a> instead');
end
