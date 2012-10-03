function res = maxlyap(ndays, disturb)
global t g_Y g_t g_grind;
i_parcheck;
if nargin == 0
   ndays = 50;
else
   ndays = i_checkstr(ndays);
end;
if nargin <= 1
   disturb = 1E-4;
else
   disturb = i_checkstr(disturb);
end;
step = 0.01;  % stepsize h in Sprott's description
N0 = i_initvar;
oldtstep = g_grind.tstep;
g_grind.tstep = ndays / step;
try
   i_ru(g_grind.odefile, t, ndays, N0, 0);
   Y = g_Y;
   times = g_t;
   xb0 = N0 + disturb;
   d0 = sqrt(size(g_Y, 2) * disturb^2);
   if ~g_grind.solver.isdiffer
      g_grind.tstep = 2;
   end;
   suml = 0;
   lambdas=zeros(1,length(times)-1);
   for i = 1:length(times) - 1
      xa1 = Y(i + 1, :)';                  % the orbit a
      i_ru(g_grind.odefile, times(i),step, xb0, 0);
      xb1 = g_Y(size(g_Y, 1), :)';         
      d1 = sqrt(sum((xb1 - xa1).^2));      % distance between both runs
      suml = suml + ln(abs(d1 / d0));      % sum of lyapunovs
      lambdas(i) = suml / i;               % running average lyapunov
      xb0 = xa1 + d0 .* (xb1 - xa1) ./ d1; % next initial condition
      d0 = sqrt(sum((xb0 - xa1).^2));      % initial distance (should not change)
   end;
   lambdas = lambdas ./ step;
   g_grind.tstep = oldtstep;
catch err
%   err=lasterror;
   g_grind.tstep = oldtstep;
   rethrow(err);
end;
i_makefig('lyap1');
%lambdas=lambdas./step;
plot(times(1:length(times) - 1), lambdas);
i_plotdefaults;
xlabel('time');
ylabel('\lambda_{max}');
res=lambdas(length(lambdas));
fprintf('Maximum lambda = %g\n', res);
