function [result, m, meanperiod, index, speciesno] = nmaxima(nsteps, silent)
global g_Y g_t t g_grind;
i_parcheck;
RelTol = 0.01;
if nargin  < 2
   silent = 0;
else
   silent = i_checkstr(silent);
end;
if nargin < 1
   nsteps = g_grind.ndays;
else
   nsteps = i_checkstr(nsteps);
end;
oldndays = g_grind.ndays;
if nsteps < 1E-30
   error('GRIND:nmaxima:TooFewTimeSteps','Number of time steps must be larger than 0');
end;
try
   g_grind.ndays = nsteps;
   N0 = i_initvar;
   if i_settingschanged(N0)
      i_ru(g_grind.odefile, t, g_grind.ndays, N0, 1);
   end;
   g_grind.ndays = oldndays;
catch err
%   err=lasterror;
   g_grind.ndays = oldndays;
   rethrow(err);
end;
if g_grind.solver.isdiffer
   tt = (0:1:nsteps)';
   Y = interp1(g_t, g_Y, tt);
else
   tt = (0:0.5:nsteps)';
   Y = interp1(g_t, g_Y, tt);
end;
iX = find(max(mean(Y)) == mean(Y));
speciesno = iX;
%to find the maxima of this species:
index = find( diff( sign( diff([0; Y(:, iX); 0]) ) ) < 0 );
if length(index) > 2
   index = index(1:length(index) - 1);
end;
if (length(index) > 2) && (index(1) == 1)
   index = index(2:length(index));
end;
meanperiod = -1;
maxima =  Y(index, iX);
index2 =  diff( sign( diff([0; Y(:, iX); 0]) ) ) > 0;
minima =  Y(index2, iX);
times_of_maxima = tt(index);
m = sort(Y(index, iX));
meanY = mean(Y(:, iX));
if ~isempty(maxima) && ~isempty(minima) && ((mean(maxima) - mean(minima)) / meanY > 0.001)
   if max(m) - 0.0001 > min(m)
      m1 = floor(m / (max(m) - min(g_Y(:, iX))) / RelTol) * RelTol;
      number_of_maxima = sum(diff(m1) > RelTol * 2) + 1;
   else
      number_of_maxima = 1;
   end;
   if (number_of_maxima == 1) && (length(times_of_maxima) > 3)
      meanperiod = mean(diff(times_of_maxima(1:length(times_of_maxima) - 1)));
   end;
   cycle = 1;
else
   number_of_maxima = 0;
   cycle = 0;
end
if ~silent
   if ~cycle
      disp('no cycles');
   else
      if number_of_maxima == 1
         fprintf('limitcycle with period of %g timesteps\n', meanperiod);
      elseif number_of_maxima < 10
         fprintf('(probably) complex cycle, %d maxima\n',number_of_maxima);
      else
         if number_of_maxima > 999
            s = '>';
         else
            s = ' ';
         end;
         fprintf('chaos or quasiperiodic cycle, %s%d maxima\n',s, number_of_maxima);
      end;
   end;
end;
if nargout > 0
   result = number_of_maxima;
end;




