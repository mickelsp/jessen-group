%function created by GRIND
function g_X2=curr_ode1(t,g_X1)
global b d h r;
fc=g_X1(1)^2/(h^2+g_X1(1)^2);
g_X2(1,1) = r*g_X1(1)*(1-g_X1(1))-b*g_X1(2)*fc;
g_X2(2,1) = b*g_X1(2)*fc-d*g_X1(2);
