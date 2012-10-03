%REDNOISE   Generate red noise
%   Function to generate white or red noise. To be used in differential
%   equations. Itself is a difference equation (Steele and Henderson):
%   T(t) = (1 -1/lambda) * (T(t-1) - T0) + T0 + beta*randn
%
%   In which:
%     T = variable with red noise
%     T0 = mean of red noise
%     lambda = 'period' of red noise. white: lambda=1 red: lambda>1.
%     beta = extend of noise
%     randn = normally distributed random number
%
%   Note: you can set the standard deviation of the resulting series to a certain value
%   with the following formula:(Ives et al., 2003):
%   beta=SD*sqrt(2/lambda-1/lambda^2) (in which SD is the resulting standard
%   deviation).
% 
%  
%   Usage:
%   REDNOISE(t,T0,lambda,beta) for explanation of coefficients, see above (t=time).
%   REDNOISE(t,T0,lambda,beta,iset) - If you use several independent rednoises 
%   within one set of equations, number them using an integer starting with 1 (iset).
%   The iset can be a vector with numbers, rednoise will then return a vector.
%   REDNOISE(t,T0,lambda,beta,0) - Calculate the noise once and use the same dataset 
%   afterwards (even not change when parameters change!)
%   REDNOISE(t,T0,lambda,beta,iset,deltat) - define the period deltat for which the
%   difference equation will be calculated (default=1)
%
%   See also model, externvar, dwiener

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:28 $
function T = rednoise(tnow, T0, labda, beta, isets, deltat)
global g_noise t g_grind;
if nargin < 5
   isets = 1;
end
if nargin < 6
   deltat = 1;
end
update = min(isets) > 0;
if isets == 0
   isets = 1;
end;
if length(isets) < length(labda)
   isets = isets(1):isets(1) + length(labda) - 1;
end
T = T0 * ones(length(tnow),length(isets));
labda = labda + zeros(length(isets), 1);
beta = beta + zeros(length(isets), 1);
for k = 1:length(tnow)
   if beta == 0
      g_noise{isets} = [];
      return;
   end;
   for i = 1:length(isets)
      iset = isets(i);
      n = round(g_grind.ndays / deltat) + 2;
      if (length(tnow)==1) && (isempty(g_noise) || ((tnow == t) && update) || (length(g_noise) < iset) || isempty(g_noise{iset}))
         g_noise{iset} = i_rednoise(n, labda(i), T0, beta(i));
      end
      if tnow(k) > t + round(deltat * length(g_noise{iset})) - 1
         n = 1000;
         g_noise{iset}(length(g_noise{iset}) + 1:length(g_noise{iset}) + n) = ...
            i_rednoise(n, labda(i), T0, beta(i), g_noise{iset}(length(g_noise{iset})));
      end;
      if (tnow(k) >= t)
         T(k, i) = g_noise{iset}(round((tnow(k) - t) / deltat) + 1);
      end;
   end;
end;
