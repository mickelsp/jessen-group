%MANIFOLDS   Trajectories in the direction of the eigenvectors
%   In the current phase plane the manifolds of the currently selected equilibrium
%   (see: findeq) are drawn. This is especially usefull for saddle points, 
%   as the stable manifolds of these equilibria are the separatrix.
%
%   Usage:
%   MANIFOLDS - draw all manifolds of the currently selected equilibrium
%   MANIFOLDS -stable - draw all stable manifolds (=separatrix).
%   MANIFOLDS -unstable - draw all unstable manifolds
%
%   See also findeq, perturb, ru, backw

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:27 $
function manifolds(stable)
global g_grind;
persize = 0.05;
if nargin == 0
   stable = 2;
else
   if strcmpi(stable(1:2), '-s')
      stable = 1;
   elseif strcmpi(stable(1:2), '-u')
      stable = 0;
   else
      stable = 2;
   end;
end;
[N1] = findeq(0);
[Jacobian, eigenvalues,  eigenvect] = i_eigen(1);
domeigen=find(eigenvalues == max(eigenvalues));
[isstable, issaddle]  = i_stability(eigenvalues, g_grind.solver.isdiffer);
if issaddle
   for i = 1:length(eigenvalues)
      stab = i_stability(eigenvalues(i), g_grind.solver.isdiffer);
      if stab && (stable~=0)
         i_keep(N1 + eigenvect(:, i) *  persize);
         backw;
         i_keep(N1 - eigenvect(:, i) *  persize);
         backw;
      elseif ~stab && (stable~=1)
         i_keep(N1 + eigenvect(:, i) *  persize);
         ru;
         i_keep(N1 - eigenvect(:, i) *  persize);
         ru;
      end;
   end;
else
   if isstable
      domeigen=find(eigenvalues == min(eigenvalues));
      i_keep(N1 + real(eigenvect(:, domeigen)) *  persize*10);
      backw;
      i_keep(N1 - real(eigenvect(:, domeigen)) *  persize*10);
      backw;
      disp('Can only show manifold of the lowest eigenvalue')
   else
      i_keep(N1 + real(eigenvect(:, domeigen)) *  persize);
      ru;
      i_keep(N1 - real(eigenvect(:, domeigen)) *  persize);
      ru;
      disp('Can only show manifold of the dominant eigenvalue')
   end;
end;



