%LYAPSPECT   Calculate Lyapunov exponent spectrum
%   Calculate the spectrum of Lyapunov exponents (lambda) and the Lyapunov dimension 
%   (=Kaplan-Yorke dimension). These parameters expresses the sensitivity to initial conditions. 
%   If there is one (clearly) positive Lyapunov exponent the model is chaotic. More positive 
%   Lyapunov exponents means "hyperchaos". The algorithm is based on: 
%   A. Wolf, J. B. Swift, H. L. Swinney, and J. A. Vastano, "Determining Lyapunov Exponents 
%   from a Time Series," Physica D, Vol. 16, pp. 285-317, 1985.
%   The Jacobian is approximated numerically if the user has not entered the equations 
%   (see enterjac)
%
%
%   Usage:
%   LYAPSPECT - calculate lyapunov spectrum with g_grind.ndays days. Default step = 0.1, StepJac=0.3
%   LYAPSPECT N STEP STEPJAC - calculate N days. The step for QR algorithm=0.1. The Jacobian is updated
%   each STEPJAC*STEP time steps. 
%   [LYAPEXP,LYAPDIM,LYAPOVS]=LYAPSPECT(N,STEP,STEPJAC) = return the lyapunov exponents (LYAPEXP) 
%   the lyapunov dimension (LYAPDIM) and all lyapunov values (LYAPOVS).
%
%   See also lyapunov, enterjac, lorenzmap, takens

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:27 $
function [lambda1, LD, lambdas,lyaps,times] = lyapspect(ndays,step,stepjac)
global g_grind t;
%tic;
if nargin < 1
   ndays = g_grind.ndays;
else
   ndays = i_checkstr(ndays);
end; 
if nargin<2
   step=1;
else
   step = i_checkstr(step);
end;
if nargin<3
   stepjac = 0.3*step;
else
   step = i_checkstr(step);
end;
N0 = i_initvar;
if g_grind.solver.nonautonomous
   if strcmp(questdlg('Might not find correct lyapunov exponents in nonautonomous equation, continue?','Error','Yes','No','No'),'No')
      error('GRIND:lypaspect:NonAutonomous','Make equation autonomous by adding tau''=1 and replacing all t by tau');
   end;
end;
d = g_grind.statevars.dim;
L = d + d^2;
A = eye(d);
IC = [N0(:); A(:)];
T = t - step;
g_grind.lyapspect.L=d^2;
g_grind.lyapspect.dtJac=stepjac*step;
g_grind.lyapspect.d=d;
g_grind.lyapspect.num=isempty(g_grind.Jacobian);
g_grind.lyapspect.J=i_calcjac(g_grind.lyapspect.num,1,N0);
g_grind.lyapspect.tJac=t+g_grind.lyapspect.dtJac;
sumR=zeros(d,1);
lyaps=zeros(ceil(ndays / step),d);
times=zeros(1,ceil(ndays / step)-1);
lambdas=zeros(ceil(ndays / step)-1,d);
for i = 1:ceil(ndays / step)
   T = T + step;
   [ts,X]=feval(str2func(g_grind.solver.name), str2func('i_runjac'),[T,T+0.5*step,T+step],IC');
   IC = X(size(X, 1), :);
   P = reshape(IC(d + 1:L), d, d);
   %QR decomposition
   [A, R] = qr(P);
   IC(d + 1:L) = A(:);
   R = abs(diag(R));
   for j = 1:size(A, 1)
      if R(j) < 1E-40
         R(j) = 0;
      else
         R(j) = log(R(j));
      end;
   end;
   lyaps(i, :) = R';
   sumR=sumR+R;
   if T>t
      lambdas(i-1,:)= sumR'/T;
      times(i-1)=T;
   end;
end;
figure;
plot(times,lambdas);
i_plotdefaults;
xlabel('time');
ylabel('\lambda''s');
lyaps=lyaps./T;
lambda = sum(lyaps);
Lambda = sort(lambda,'descend');
%To calculate the Lyapunov dimension (or Kaplan-Yorke dimension)
LESum = Lambda(1); 			
LD = 0;
if (d > 1 && Lambda(1) > 0)
   for N = 1:d - 1
      if Lambda(N + 1) ~= 0
         LD = N + LESum / abs(Lambda(N + 1));
         LESum = LESum + Lambda(N + 1);
         if LESum < 0
            break;
         end
      end
   end
end;
for i = 1:d
   fprintf('lambda(%g) = %g\n', i, lambda(i));
end;
fprintf('Lyapunov dimension (or Kaplan-Yorke dimension) = %g\n', LD);
if nargout > 0
   lambda1 = lambda;
end;
g_grind.lyapspect = [];
%toc;
