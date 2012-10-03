function dist = i_totgrowth2(x)
% version for difference equations
global g_grind;
x0=x;
for i=1:g_grind.solver.iters
   x0=feval(g_grind.odefile, 1, x0);
end;
Nres = x0 - x;
dist = sum(Nres.^2);

