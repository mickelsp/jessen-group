%function created by GRIND
function g_X2=curr_ode4(t,g_X1)
global K r;
g_X2(1,1) = g_X1(1)*r*(1-g_X1(1)/K);
