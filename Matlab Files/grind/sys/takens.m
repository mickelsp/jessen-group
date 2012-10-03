%TAKENS   Create a Takens plot 
%   Plot of N(t) versus N(t+timestep). TAKENS analyses the results of the 
%   last run. (if there is no last run or if parameters have changed, it calls RU). 
%   Use ru if you want to update the last run.
%
%   Usage:
%   TAKENS uses a default time lag  of 1 time step.
%   TAKENS X uses a time lag of X time steps.
% 
%
%
%   See also lorenzmap, lyapunov

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:29 $
function takens(timestep)
global g_Y g_t t g_grind;
i_parcheck;
if nargin == 0
   timestep = 1;
else
   timestep = i_checkstr(timestep);
end;
N0 = i_initvar;
if i_settingschanged(N0)
   i_ru(g_grind.odefile, t, g_grind.ndays, N0, 1);
end;
lag=interp1(g_t,g_Y,g_t+timestep);
i_makefig('takens');
oldhold = ishold;
hold on;
leg = cell(g_grind.statevars.dim, 1);
Hs = zeros(g_grind.statevars.dim, 1);
for iX = 1:size(g_Y, 2)
   hp = plot(g_Y(:, iX), lag(:, iX), '.');
   set(hp, 'Color', g_grind.pen.color2);
   Hs(iX)=hp;
   leg{iX}=i_statevars_names(iX);
   nextpen;
end;
legend(Hs,leg);
xlabel('X_t');
ylabel(['X_{t+' num2str(timestep) ')}']);
if ~oldhold
   hold off;
end;





