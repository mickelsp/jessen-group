%POINCAREMAP   Construct a Poincare map
%   All trajectories that cross a certain surface in the state variable space
%   are mapped to itself, i.e. the subsequent crossings are mapped.
%   POINCAREMAP analyses the results of the last run. (if there is
%   no last run, it calls time). Use RU if you want to
%   update the last run.
%
%   Usage:
%   POINCAREMAP VAR1 VAR2 VALUEVAR2 - makes a map of VAR1 on the plane VAR2=VALUEVAR2. 
%   Only increasing trajectories are analysed.
%   POINCAREMAP VAR1 VAR2 VALUEVAR2 0 - Analyses only decreasing trajectories.
%  
%   See also POINCARESECT, LORENZMAP

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:28 $
function poincaremap(avar1, avar, avalue, increasing)
%global g_Y;
i_parcheck;
if nargin < 3
   errordlg('Not enough arguments, Usage: poincaremap avar1 avar avalue');
   error('GRIND:poincaremap:ArgError','Not enough arguments, Usage: poincaremap avar1 avar avalue');
end;
if nargin < 4
   increasing = 1;
else
   increasing = i_checkstr(increasing);
end;
avalue = i_checkstr(avalue);
us.avar1 = avar1;
us.avar = avar;
us.avalue = avalue;
us.increasing = increasing;
poincar = i_poincare(avar, avalue, increasing);
ivar = i_varno(avar1);
H = i_makefig('poinmap');
set(H,'name',['Poincare map through plane ' avar ' = ' num2str(avalue)])
set(H, 'userdata', us);
m = i_lagmap(poincar);
plot(poincar(:, ivar), m(:, ivar), 'k.');
xlabel(i_disptext([avar1 '(t)']));
ylabel(i_disptext([avar1 '(t+1)']));


