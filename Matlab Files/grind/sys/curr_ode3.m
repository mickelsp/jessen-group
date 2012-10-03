%function created by GRIND
function g_X2=curr_ode3(t,g_X1)
global c h K p r;
g_X2(1,1) = r*g_X1(1)*(1-g_X1(1)/K)-c*g_X1(1)^p/(g_X1(1)^p+h);
