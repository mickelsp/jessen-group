%PERTURB   Perturb saddle point
%   perturb a saddle point in the direction of the attracting eigenvector
%   use BACKW after perturbation of the equilibrium to approximate one part of
%   the separatrix
%
%   To find the whole separatrix:
%   perturb 1
%   backw 30
%   perturb -1
%   backw 30
%
%   Usage:
%   PERTURB - Perturb the equilibrium in positive direction with a size of
%   0.05
%   PERTURB -1 SIZE - direction can be 1 or -1. SIZE should be a small
%   number.
%
%   See also backw, findeq

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:28 $
function perturb(dir, persize)
global g_grind;
i_parcheck;
if nargin == 0
   dir = 1;
else
   dir =i_checkstr(dir);
end;
if nargin < 2
   persize = 0.005;
else
   persize = i_checkstr(persize);
end;
[N1] = findeq(0);
[Jacobian, eigenvalues, eigenvect] = i_eigen(1);
if abs(dir)>10
   vec=abs(dir)-10;
   dir=dir/abs(dir);
else
   if i_stability(eigenvalues(1),g_grind.solver.isdiffer)
     vec=1;
   elseif i_stability(eigenvalues(2),g_grind.solver.isdiffer)
      vec=2;
   else
      vec=1;
   end;
end;
if abs(dir) == 2
   dir = dir / 2;
   if vec == 1
      vec = 2;
   else
      vec = 1;
   end;
end;
len = length(real(eigenvect(:, vec)));
if len < 0.1
   len = 0.1;
end;
N0 = N1 + dir * real(eigenvect(:, vec)) / len * persize;
i_keep(N0);

