%i_EIGEN = calc Jacobian and eigenvalues of the current initial condition
%
function [Jacobian, eigenval, eigenvect] = i_eigen(donumerical,niters,N0)
global g_grind;
if nargin<1
   donumerical=isempty(g_grind.Jacobian);
end;
if nargin < 3
   N0 = i_initvar;
end;
if nargin < 2
   niters = g_grind.solver.iters;
end;
if isempty(g_grind.Jacobian)
   donumerical = 1;
end;
Jacobian = i_calcjac(donumerical, niters,N0);
[eigenvect, eigenval] = eig(Jacobian);
eigenval = diag(eigenval);
%eigenval=eig(Jacobian);
