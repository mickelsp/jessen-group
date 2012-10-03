function H = i_figno(afig)
global g_grind;
if ~isfield(g_grind, 'timevars')
   ntime = 8;
else
   ntime = max(8, length(g_grind.timevars));
end;
switch afig
 case 'time'
   H = 1;
 case 'phase2'
   H = 1 + ntime;
 case 'phase3'
   H = 2 + ntime;
 case 'phase1'
   H = 3 + ntime;
 case 'funplot'
   H = 4 + ntime;
 case 'dirfield'
   H = 5 + ntime;
 case 'eigen'
   H = 6 + ntime;
 case 'paranal'
   H = 7 + ntime;
 case 'poinsec'
   H = 10 + ntime;
 case 'potent1'
   H = 11 + ntime;
 case 'potent2'
   H = 12 + ntime;
 case 'potent3'
   H = 13 + ntime;
 case 'poinmap'
   H = 14 + ntime;
 case 'lyap1'
   H = 15 + ntime;
 case 'lyap2'
   H = 16 + ntime;
 case 'itermap'
   H = 17 + ntime;
 case 'takens'
   H = 18 + ntime;
 case 'lorenzmap'
   H = 19 + ntime;
 case 'combfig'
   H = 20 + ntime;
 case 'autocorr'
   H = 21 + ntime;
 case 'torus'
   H = 22 + ntime;
 case 'attrbasin'
   H = 23 + ntime;
 case 'trdet'
   H = 24 + ntime;
 case 'dialog'
   H = 25 + ntime;
 case 'setmat'
   H = 26 + ntime;
 case 'mcarlo'
   H = 27 + ntime; %4 figs
 case  'conteq2d'
   H = 31 + ntime;
 case  'conteq'
   H = 32 + ntime; 
case  'growths'
   % for each state variable one growth plot
   H = 33 + ntime;
 case 'obspred'
   H = 34 + ntime+ nstatevars;
 case 'vectplot'
   H = 35 + ntime + 2 * nstatevars;
 case 'varcontour' %obsolete, synonym of vectplot
   H = 35 + ntime + 2 * nstatevars;
 case 'uncertain'
   H = 35 + ntime + 3 * nstatevars;
 case 'paranal2d'
   H = 35 + ntime + 4 * nstatevars;
 case  'viewcells';
   H = 35 + ntime + 4 * nstatevars + 1 * nvectors;
 case 'maxno'
   H = 35 + ntime + 4 * nstatevars + 2 * nvectors + 1;
 otherwise
   error('GRIND:figno:UnknownType','Unknown figure: %s', afig);
end;
function n = nstatevars
global g_grind
if ~isfield(g_grind, 'statevars')
   n = 0;
else
   n = g_grind.statevars.dim;
end;
function n = nvectors
global g_grind;
if ~isfield(g_grind, 'statevars')
   n = 0;
elseif ~isempty(g_grind.statevars) && g_grind.statevars.vector
   n = length(g_grind.statevars.dims);
else
   n = 0;
end;
