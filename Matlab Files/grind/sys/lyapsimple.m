%LYAPSIMPLE - Simple Elner algoritm for lyapunov exponent
%   Use simtime to get a time series with a small fixed time step
%
%   Usage:
%   lyapsimple - gives the lyapunov

function LE=lyapsimple(timestep)
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
LE = 0;
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
   LE = LE + log(maxU);
   U = U / maxU;
end;
LE = LE / (tt(end) - tt(1));
disp('LE = ');
disp(LE);
