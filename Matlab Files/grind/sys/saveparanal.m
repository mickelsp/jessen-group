function saveparanal(matfile)
global g_paranal g_grind; %#ok
if isfield(g_grind,'paranal2d')
   answer=g_grind.paranal2d; %#ok
elseif isfield(g_grind,'paranal')
   answer=g_grind.paranal; %#ok
else
   error('GRIND:saveparanal:NoRun','Cannot save paranal, not yet run');
end;
save(matfile,'g_paranal','answer','g_grind');
