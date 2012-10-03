function dist = i_totgrowth(x)
global g_grind;
Nres = feval(g_grind.odefile, 1, x);
dist = sum((Nres).^2);