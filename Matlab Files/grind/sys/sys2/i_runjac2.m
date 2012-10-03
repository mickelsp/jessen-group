function res=i_runjac2(t,g_A)
global g_grind;
L=length(g_A);
d=sqrt(L);
%linear interpolation of the Jacobian matrix
J=g_grind.lyapspect.prevJ+(g_grind.lyapspect.J-g_grind.lyapspect.prevJ)...
    *t/g_grind.lyapspect.t;
res=reshape(J*reshape(g_A,d,d),L,1);

