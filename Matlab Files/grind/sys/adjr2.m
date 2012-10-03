%ADJR2   Calculate mean adjusted R2
%   The adjusted R2 is a measure of the fit:
%   it is defined as:
%   AdjR2 = 1 - (UnexplainedSS * (n - 1)) / (TotalSS * (n - NEstimPars - 1));
%
%   Usage:
%   ADJR2(NPARS) - Calculate adjusted R2. NPARS is the number of parameters used 
%   for fitting.
%
%   See also optimpars

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:25 $
function meanAdjR2 = adjr2(npars)
global g_data g_grind t g_Y g_t;
if (nargin<1)||isempty(g_data.pred)
   if isempty(g_data)
      error('GRIND:adjr2:NoData','There are no data for calculating adjusted R2, use SETDATA to add data');
   end;
   if isfield(g_data,'pars')&&~isempty(g_data.pars)
      npars=length(g_data.pars);
   else
      npars=1;
   end;
   i_ru(g_grind.odefile, t, g_data.t(length(g_data.t),1), i_initvar,0);
   for i = 1 :g_grind.statevars.dim
     g_data.pred(:, i) = interp1(g_t, g_Y(:, i), g_data.t);
   end;
end;
if isempty(g_data.obs)
   error('GRIND:adjr2:NoData','no data entered');
end;
g_t = [];
g_Y = [];
par;
val;
disp(' ');
totAdjR2 = 0;
% ifig = i_figno('obspred');
nn = 0;
for i = 1:g_grind.statevars.dim
   if ~isnan(g_data.minobs(i))
      i_makefig('obspred', i);
      nn = nn + 1;
      plot(g_data.obs(:,i),g_data.obs(:,i),'k',g_data.obs(:,i),g_data.pred(:,i), 'o');
      xlabel(i_disptext(['Observed ' i_statevars_names(i)]));
      ylabel(i_disptext(['Predicted ' i_statevars_names(i)]));
      AdjR2 = GetAdjR2(g_data.pred(:, i), g_data.obs(:, i), npars);
      totAdjR2 = totAdjR2 + AdjR2;
      title(['Adjusted R^{2} = ' num2str(AdjR2)]);
   end;
end;
disp(['Mean adjusted R2  : ' num2str(totAdjR2 / nn)]);
if nargout > 1
   meanAdjR2 = totAdjR2;
end;

function [AdjR2, TotalSS, UnexplainedSS, R2] = GetAdjR2(Pred, Obs, NEstimPars)
TotalSum = 0;
n = 0;
nObs = length(Obs);
for i = 1:nObs
   if ~isnan(Obs(i))
      n = n + 1;
      TotalSum = TotalSum + Obs(i);
   end;
end;
if n > 0
   Mean = TotalSum / n;
   TotalSS = 0;
   UnexplainedSS = 0;
   for i = 1:nObs
      if ~isnan(Obs(i))
         TotalSS = TotalSS + (Obs(i) - Mean)^2;
         UnexplainedSS = UnexplainedSS + (Obs(i) - Pred(i))^2;
      end;
   end;
   % Calculate the Adjusted R2
   if (n - NEstimPars > 1) && (TotalSS > 0)
      AdjR2 = 1 - (UnexplainedSS * (n - 1)) / (TotalSS * (n - NEstimPars - 1));
      R2=1-UnexplainedSS/TotalSS;
   elseif (TotalSS > 0) && (NEstimPars == 0)
      AdjR2 = 1 - UnexplainedSS / TotalSS;
      R2=AdjR2;
   else
      AdjR2 = 0;
      R2=0;
   end;
else
   AdjR2 = NaN;
end;
