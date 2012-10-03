%solve an odefunct as a difference equation
%
function [T, Y] = i_differ(odefunct, TSPAN, Y0, OPTIONS)
global g_grind g_Y g_t;
if isempty(g_grind)
   iters = 1;
else
   iters = g_grind.solver.iters;
end;
if length(TSPAN) == 2
   T = transpose(TSPAN(1):TSPAN(2));
else
   T=TSPAN;
   if size(T,1)<size(T,2)
      T = transpose(T);
   end;
end;
if size(Y0,2)>size(Y0,1)
   Y0=Y0';
end;
T0 = T(1);
ndays = T(size(T, 1)) - T0;
daily=ndays == size(T, 1) - 1;
if daily
   Y = zeros(size(T, 1), size(Y0, 1));
   N0 = Y0;
   Y(1, :) = transpose(N0);
   for i = 2:ndays + 1
      for k = 1:iters
         N0 = feval(odefunct, T0 + i, N0);
         if max(isnan(N0)) && g_grind.solver.backwards
%            global g_Y g_t;
            g_Y = Y(1:i - 1, :); %#ok
            g_t = T(1:i - 1); %#ok
            error('GRIND:backw:differ','Backwards: error solving difference equation backwards');
         end;
      end;
      Y(i, :) = transpose(N0);
   end;
else
   lenT = size(T, 1);
   Y = zeros(lenT, size(Y0, 1));
   Y(1, :) = transpose(Y0);
   N0 = Y0;
   for iT = 2:lenT
      NN = round(T(iT) - T0);
      for j = 1:NN
         for k = 1:iters
            N0 = feval(odefunct, T0 + j, N0);
         end;
      end;
      T0 = T0 + NN;
      Y(iT, :) = transpose(N0);
   end;
end;




