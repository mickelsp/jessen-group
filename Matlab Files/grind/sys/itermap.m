%ITERMAP   Iteration map for 1-D difference equation
%   Plot x(t) versus x(t+n) for a 1-D difference equation (default
%   n=1).
%
%   Usage: 
%   ITERMAP - plot x(t) versus x(t+1). This is simply a plot of
%   the right-hand-sides of the difference equation (with the
%   diagonal).
%   ITERMAP N - plot x(t) versus x(t+N) of the state variable
%   on the x-axis (see ax). N can also be
%   negative for backwards simulations.
%   ITERMAP [N1 N2 .. Nn] - you can plot more than one itermaps at a time
%   by replacing the scalar N by a vector with the N's.
%   ITERMAP N VAR LIM - plot VAR(t) versus VAR(t+N) with a range
%   of LIM.
%   
%   See also ru, null, ax

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:27 $
function itermap(n, avar, alim)
global g_grind t;
i_parcheck;
if nargin == 0
   n = g_grind.solver.iters;
else
   n = i_checkstr(n);
end;
if nargin<=1
   avar=g_grind.xaxis.var;
end;
if nargin<=2
   alim=g_grind.xaxis.lim;
else
   alim=i_checkstr(alim);
end;
if ~g_grind.solver.isdiffer
   err = i_warningdlg('GRIND:itermap:nodifference','This function is designed for difference equations, please use <a href="matlab:plotdiff">plotdiff</a> instead');
else
   err = -1;
end
if length(n)>1
   [H,new] = i_makefig('phase1');
   set(H, 'Name', 'Iteration map');
   if new
      plot(alim,alim,'k');
   end;
   leg=cell(length(n)+1,1);
   leg{1}='y=x';
   for i=1:length(n)
      itermap(n(i),avar,alim);
      leg{i+1}=i_disptext([char(avar) '_{t+' num2str(n(i)) '}']);
   end;
   if new
      legend(leg,2);
   end;
   yLabel([i_disptext(char(avar)) '_{t+n}']);
   return;
end;
N0 = i_initvar;
iX = i_getno(avar);
if ~iX.isvar
   errordlg('Need to have a state variable on the x-axis');
   error('GRIND:itermap:NoStatevars','Need to have a state variable on the x-axis');
end;
X = alim(1):(alim(2) - alim(1)) / 200:alim(2);
Y = zeros(size(X, 2), size(N0, 1));
for i = 1:size(X, 2)
   N0(iX.no) = X(i);
   N1 = N0;
   for j = 1:abs(n)
      N1 = feval(g_grind.odefile, t, N1);
   end;
   Y(i, :) = transpose(N1);
end;
[H,new] = i_makefig('phase1');
if new
set(H, 'WindowButtonDown', 'i_callb(''mdown'')');
set(H, 'WindowButtonMotionFcn', 'i_callb(''mmove'')');
end;
% if ~g_grind.version.isoctave
%   ud=get(H,'userdata');
%   ud.iters=n;
%   set(H, 'userdata',ud);
% end;
oldhold = ishold;
hold on;
if new
   plot(alim, alim, 'k');
end;
apen = g_grind.pen;
apen.i = n - 1;
apen = nextpen(apen);
if n>0
   plot(X, Y, 'Color', apen.color2);
elseif n<0
   plot(Y, X, 'Color', apen.color2);
else
   plot(X,X, 'Color',apen.color2);
end;
if g_grind.statevars.dim > 1
   title(['Valid for ' i_othervars(N0, iX.no)]);
end;
xlabel(i_disptext([char(avar) '_{t}']));
ylabel(i_disptext([char(avar) '_{t+' num2str(n) '}']));
set(gca, 'XLim', alim);
if ~oldhold
   hold off
end;
if ishandle(err)
   figure(err)
end;

