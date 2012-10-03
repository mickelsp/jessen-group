%LYAPSIMPLE - Simple Elner algoritm for lyapunov exponent
%   Use simtime to get a time series with a small fixed time step
%
%   Usage:
%   lyapsimple - gives the lyapunov

function locallyap(timestep)
global g_grind t g_t g_Y;
if nargin == 0
   if g_grind.solver.isdiffer || isnan(g_grind.tstep)
      timestep = 1;
   else
      timestep =  g_grind.ndays / g_grind.tstep;
   end;
else
   timestep = i_checkstr(timestep);
end;
N0 = i_initvar;
if i_settingschanged(N0)
   i_ru(g_grind.odefile, t, g_grind.ndays, N0, 1);
end;
tt = (g_t(1):timestep:g_t(end))';
YY = interp1(g_t, g_Y, tt);
Jac = zeros(length(tt), length(N0)^2);
U = ones(size(N0));
LEs = zeros(length(tt),1);
donum = isempty(g_grind.Jacobian);
for i = 1:length(tt)
   N0 = YY(i, :)';
   J = i_calcjac(donum, 1,N0);
   Jac(i, :) = J(:)';
   if g_grind.solver.isdiffer
      U = J * U;
   else
      U = expm(timestep * J) * U;
   end;
   maxU = max(abs(U));
   LEs(i) = LEs(i) + log(maxU);
   U = U / maxU;
end;
LE = sum(LEs) / (tt(end) - tt(1));
plot(tt,LEs)
fprintf('Global Lyapunov = %g\n',LE);
fprintf('Window\t Number of LE\t  Average LE\t  STD LE\t 1st quantile LE\t Median\t 2nd quantile LE\t%%positive\n');
for window = [10, 20, 50, 100, 200, 400, 600, 800]
   if window<tt(end)
    LLE=zeros(length(tt)-window,1);
    for i=1:length(tt)-window
       LLE(i)=sum(LEs(i:i+window))/(tt(i+window)-tt(i));
    end;
    percpos=sum(LLE>0)/length(LLE)*100;
     fprintf('%4g\t %2.4f\t %2.4f\t %2.4f\t %2.4f\t %2.4f\t %2.4f %4.2f\n', ...
        window, length(LLE), mean(LLE),std(LLE),i_makepercentiles(LLE',0.25),...
        i_makepercentiles(LLE',0.50),i_makepercentiles(LLE',0.75),percpos);
    end;
   %     subplot(4, 2, index);
   %boxplot(LLE(:),'labels', int2str(window));
end;
if g_grind.ndays<1000
   warning('GRIND:locallyap:longrun', 'locallyap needs long runs, use <a href="matlab:simtime>simtime</a> to increase the run >>1000');
end;
      

