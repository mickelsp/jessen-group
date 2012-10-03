function loadparanal(f)
global g_paranal g_grind; %#ok
inif=g_grind.inifile;
load(f);
if ~strcmp(inif,g_grind.inifile)
  use(g_grind.inifile);
%  global g_paranal g_grind;
  load(f);
end;
g_grind.paranal=answer;
paranal 1