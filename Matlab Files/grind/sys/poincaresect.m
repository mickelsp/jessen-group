%POINCARESECT   Construct a Poincare section
%   All trajectories that cross a certain surface in the state variable space
%   are plotted.
%   POINCARESECT analyses the results of the last run. (if there is
%   no last run, it calls time). Use ru if you want to
%   update the last run.
%
%   Usage:
%   POINCARESECT VAR VALUEVAR - analyses the plane VAR=VALUEVAR. 
%   Only increasing trajectories are analysed.
%   POINCARESECT VAR VALUEVAR 0 - Analyses only decreasing trajectories.
%  
%   See also poincaremap, lorenzmap

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:28 $
function poincaresect(avar, avalue, increasing)
global g_grind;
i_parcheck;
if nargin < 2
   errordlg('Not enough parameters, Usage: poincaremap avar avalue');
   error('GRIND:poincaresect:ArgError','Not enough parameters, Usage: poincaremap avar avalue');
end;
if nargin < 3
   increasing = 1;
else
   increasing = i_checkstr(increasing);
end;
avalue = i_checkstr(avalue);
[poincar, ivar] = i_poincare(avar, avalue, increasing);
i = 1;
if i == ivar
   i = i + 1;
end;
j = i + 1;
if (j == ivar) && (j < g_grind.statevars.dim)
   j = j + 1;
end;
if j <= g_grind.statevars.dim
   us.avar = avar;
   us.avalue = avalue;
   us.increasing = increasing;
   H = i_makefig('poinsec');
   plot(poincar(:, i), poincar(:, j), 'k.');
   set(H,'name',['Poincare section through plane ' avar ' = ' num2str(avalue)])
   set(H, 'userdata', us);
   title(['Poincare section S: ' avar ' = ' num2str(avalue)]);
   xlabel(i_disptext(i_statevars_names(i)));
   ylabel(i_disptext(i_statevars_names(j)));
else
   error('GRIND:poincaresect:TooFewDims','Not enough dimensions to draw a 2D poincare section');
end;

