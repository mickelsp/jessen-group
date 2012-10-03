%LORENZMAP   Create a Lorenz map
%   Plot subsequent maximums of a variable.
%
%   Usage:
%   LORENZMAP VAR1 plot a Lorenz map of the variable VAR1.
%   LORENZMAP plot the variable of the x axis of the phase plane
% 
%   LORENZMAP analyses the results of the last run. (if there is
%   no last run or parameters have changed, it calls RU). Use
%   ru if you want to update the last run.
%
%   See also ru, takens, poincaremap

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:27 $
function lorenzmap(avar)
global g_Y g_grind t;
i_parcheck;
if nargin == 0
   avar = g_grind.xaxis.var;
end;
N0 = i_initvar;
if i_settingschanged(N0)
   i_ru(g_grind.odefile, t, g_grind.ndays, N0, 1);
end;
iX = i_varno(avar);
if isempty(iX)
   maxima=i_maxima(i_getoutfun(g_grind.xaxis.var),1);
else
   maxima = i_maxima(g_Y, iX);
end;
i_makefig('lorenzmap');
maxlag = i_lagmap(maxima);
hp=plot(maxima,maxlag,'.','Color',g_grind.pen.color);
set(hp, 'userdata', avar);
xlabel(i_disptext([ avar '_{max, t}']));
ylabel(i_disptext([ avar '_{max, t+1}']));
hold on;
rang = [min(maxima), max(maxima)];
plot(rang, rang, 'k');
hold off;






