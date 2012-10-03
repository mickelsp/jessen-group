function  i_setsettings(settings)
global g_grind t;
g_grind.ndays=settings(1);
t=settings(2);
g_grind.tstep=settings(3);
g_grind.solver.iters=settings(4);
g_grind.solver.backwards=settings(5);
g_grind.solver.addmode=settings(6);
if settings(7)>0
    g_grind.solver.opt.RelTol=settings(7);
end;
if settings(8)>0
   g_grind.solver.opt.AbsTol=settings(8);
end;
if settings(9)>0
   g_grind.solver.opt.MaxStep=settings(9);
end;
solverlist = {'i_differ','ode45','ode23','ode113','ode15S', 'ode23S', 'ode23T', 'ode23TB','rk4','Euler',};
g_grind.solver.name=solverlist{settings(10)};
j = 11;
for i = 1:size(g_grind.pars, 2)
   p = evalin('base', char(g_grind.pars{i}));
   s=size(p);
   jj=prod(s);
   assignin('base',char(g_grind.pars{i}),reshape(settings(j+1:j+jj),s(1),s(2)));
   j=j+jj;
end;
N0= i_initvar;
jj=numel(N0);
%jj2=numel(NP);
i_keep(settings(j+1:j+jj));
%,settings(j+jj+1:j+jj+jj2));




