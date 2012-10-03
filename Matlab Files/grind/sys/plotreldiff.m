%PLOTRELDIFF   Plot per capita growth of 1D differential equation
%   Create a plot of the per capita growth of a differential equation of a 
%   one dimensional system. On the y-axis the relative growth is plotted.
%
%   Usage:
%   PLOTRELDIFF - plots the variable of the x axis with the range of the x-axis.
%   PLOTRELDIFF VAR1 LIM - plots the state variable VAR1 with a range of LIM. 
%
%   See also plotdiff, null, phas, ax

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:28 $
function plotreldiff(avar, alim)
global g_grind;
i_parcheck;
if nargin==0
   avar=g_grind.xaxis.var;
end;
if nargin<=1
   alim=g_grind.xaxis.lim;
else
   alim=i_checkstr(alim);
end;
%if size(g_grind.statevars,2)>1
%   errordlg('Warning: This command is designed for 1D differential equations');
%end;
iX = i_varno(avar);
if isempty(iX)
   error('GRIND:plotreldiff:NoStatevar','Can only create plot if state variable is on the axis');
end;
[H,new] = i_makefig('funplot');
set(H, 'Name','Plot of per capita growth')
if new
set(H, 'WindowButtonDown', 'i_callb(''mdown'')');
set(H, 'WindowButtonMotionFcn', 'i_callb(''mmove'')');
end;
N0 = i_initvar;
if alim(1)==0 
   alim(1)=(alim(2) - alim(1)) / 1000;
end;
X = (alim(1):(alim(2) - alim(1)) / 1000:alim(2));
Y = zeros(1, size(X, 2));
for i = 1:size(X, 2)
   N0(iX) = X(i);
   yy = feval(g_grind.odefile, 1, N0);
   Y(i) = yy(iX)/X(i);
end;
oldhold = ishold;
hold on;
plot(X, Y);
xlabel(i_disptext(avar));
ylabel([i_disptext(avar) '''/' i_disptext(avar)]);
plot(X, zeros(1, size(X, 2)), 'k');
if g_grind.statevars.dim > 1
   title(i_disptext(['Valid for ' i_othervars(N0, iX)]));
end;
set(gca, 'XLim', alim);
if ~oldhold
   hold off;
end;
if g_grind.solver.isdiffer
   i_warningdlg('GRIND:plotreldiff:differenceeq','This function is designed for differential equations, please use <a href="itermap">itermap</a> instead');
end
