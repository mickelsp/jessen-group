function changed = i_settingschanged(N0, ndays)
global g_grind g_Y g_t;
if nargin == 0
   settings = i_getsettings;
else
   settings = i_getsettings(N0);
end;
if nargin == 2
   settings(1) = ndays;
end;
if isempty(g_grind.lastsettings)
   changed = 1;
elseif g_grind.solver.addmode
   % in addmode always run, update the run if the initial conditions change
   nperm=0;
   if isfield(g_grind,'permanent')
    for i=1:size(g_grind.permanent)
     p = evalin('base', char(g_grind.permanent{i}.name));
     nperm=nperm+numel(p);
    end;
   end;
   if (nargin>0) && isdifferent(N0,g_grind.lastsettings(end-length(N0)+1-nperm:end-nperm)) || (size(g_Y, 1) < 3) || (size(g_Y, 1)~=length(g_t));
      g_Y=[];
      g_t=[];
   end;
   changed=1;
elseif isempty(g_t) || ((nargin == 2) && (abs(g_t(length(g_t)) - settings(1) + settings(2)) > 1))
   changed = 1;
else
   changed = isdifferent(settings,g_grind.lastsettings) || (size(g_Y, 1) < 3);
end;

function res = isdifferent(A, B)
res=~min(size(A)==size(B));
if isempty(B)&&isempty(A)
   res = 0;
elseif ~res
   compar = A == B;
   res = ~min(compar);
   if res
      res = ~min(compar + isnan(A) .* isnan(B));
   end;
end;
