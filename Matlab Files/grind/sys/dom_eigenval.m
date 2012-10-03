%DOM_EIGENVAL - dominant (maximum) eigenvalue
%
%  Usage:
%  Use DOM_EIGENVAL in OUT or AX
%
function res = dom_eigenval
global g_Y g_grind;
res=ones(size(g_Y,1),1);
for i=1:size(g_Y,1)
   N0=g_Y(i,:)';
   [T1,V1]=i_eigen(isempty(g_grind.Jacobian),g_grind.solver.iters,N0);
   res(i)=max(real(V1));
end;

   