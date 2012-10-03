function [settings] = i_getsettings(N0)
global g_grind t;
maxn = 128;
settings = zeros(maxn,1);
settings(1)=g_grind.ndays;
settings(2)=t;
settings(3)=g_grind.tstep;
settings(4)=g_grind.solver.iters;
settings(5)=g_grind.solver.backwards;
settings(6)=g_grind.solver.addmode;
if ~isempty(g_grind.solver.opt.RelTol)
   settings(7)=g_grind.solver.opt.RelTol;
end;
if ~isempty(g_grind.solver.opt.AbsTol)
   settings(8)=g_grind.solver.opt.AbsTol;
end;
if ~isempty(g_grind.solver.opt.MaxStep)
   settings(9)=g_grind.solver.opt.MaxStep;
end;
solverlist = {'i_differ','ode45','ode23','ode113','ode15S', 'ode23S', 'ode23T', 'ode23TB','rk4','Euler',};
for i = 1:length(solverlist)
   if strcmpi(solverlist{i}, g_grind.solver.name)
      settings(10) = i;
      break;
   end;
end;
j = 11;
for i = 1:length(g_grind.pars)
   p = evalin('base', char(g_grind.pars{i}));
   jj=numel(p);
   j2=j+jj;
   while j2 + 1 > maxn
      addn=j2 + 1;
      maxn = maxn + addn;
      settings = [settings; zeros(addn,1)];
   end;
   settings(j+1:j2)=p(1:jj);
   j=j2;
end;
if nargin == 0
   N0 = i_initvar;
end;
settings=[settings(1:j);N0];
j=j+length(N0);
if isfield(g_grind,'permanent')
    for i=1:size(g_grind.permanent)
     p = evalin('base', char(g_grind.permanent{i}.name));
     jj=numel(p);
     j2=j+jj;
     while j2 + 1 > maxn
        addn=j2 + 1;
        maxn = maxn + addn;
        settings = [settings; zeros(addn,1)];
     end;
     settings(j+1:j2)=p(1:jj);
     j=j2;
    end;
end;



