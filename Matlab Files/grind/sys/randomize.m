%Randomize the random generator
%shorthand for:
%   rand('state',sum(100*clock));
%global g_grind;
%if ~isempty(g_grind)&&~isnan(g_grind.randstate)
%   rand('state',g_grind.randstate);
%else
%   rand('state',sum(100*clock));
%end;
 rng('default')