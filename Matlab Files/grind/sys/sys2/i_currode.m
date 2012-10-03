function res=i_currode(at,ay,yp)
global g_grind;
newyp=feval(str2func(g_grind.odefile),at,ay);
res=yp-newyp;