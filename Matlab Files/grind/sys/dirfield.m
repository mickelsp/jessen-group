%DIRFIELD   Create direction field of a simple differential equation
%   Create a direction field for studying single first order
%   differential equations. You can click in the figure to run 
%   from that point.
%
%   Usage:
%   DIRFIELD - Create a direction field for 100 days and for the
%   variable which was selected for the X axis (see also ax).
%   DIRFIELD VAR1 LIM - Create a direction field for VAR1 with a
%   range of LIM.
%   
%   Examples:
%   DIRFIELD X [0 10]
%
%   See also ru, time, null

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:26 $
function dirfield(avar, alim)
global g_grind t;
i_parcheck;
npoints = 35;
if g_grind.ndays > 100
   ndays = 100;
else
   ndays = g_grind.ndays;
end
if nargin == 0
   avar = g_grind.xaxis.var;
end;
if nargin < 2
   alim = g_grind.xaxis.lim;
end;
iY = i_getno(avar);
if isempty(iY.no)
   errordlg('Cannot create direction field if there are no state variables on the ''x''axis.');
   error('GRIND:dirfield:NoStateVariable','Cannot create direction field if there are no state variables on the ''x''axis.');
end;
Xaxis = [t, t + ndays];
Yaxis = alim;
N = i_initvar;
Vect = zeros(npoints, npoints, 2);
X = zeros(npoints, npoints);
Y = zeros(npoints, npoints);
isdiffer = g_grind.solver.isdiffer;
for y1 = 1:npoints
   N(iY.no) = (y1 - 1) * (Yaxis(2) - Yaxis(1)) / npoints + Yaxis(1);
   for t1 = 1:npoints
      tnow = (t1 - 1) * ndays / npoints + t;
      Nres = feval(g_grind.odefile, tnow, N);
      X(t1, y1) = tnow;
      Y(t1, y1) = N(iY.no);
      if isdiffer
         Nres = Nres - N;
      end
      Vect(t1, y1, 2) = Nres(iY.no);
      Vect(t1, y1, 1) = 1;
      %g_grind.ndays/(npoints*2);
   end;
end;
[H,new] = i_makefig('dirfield');
if new
   set(H, 'WindowButtonDown', ['i_callb(''mdown2'',''' avar ''',' num2str(ndays) ')' ]);
   set(H, 'WindowButtonMotionFcn', 'i_callb(''mmove'')');
%plotedit off;
end;
set(H,'Name','Direction field');
oldhold = ishold;
hold on;
h = quiver(X, Y, Vect(:, :, 1), Vect(:, :, 2));
set(h, 'Color', [0.5 0.5 0.5]);
set(gca, 'XLim', Xaxis);
set(gca, 'YLim', alim);
xlabel('t');
ylabel(i_disptext(avar));
if ~oldhold
   hold off;
end;

