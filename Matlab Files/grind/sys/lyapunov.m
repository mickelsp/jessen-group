%LYAPUNOV   The maximal Lyapunov exponent
%   Calculate the maximal Lyapunov exponent (lambda) by running a model twice 
%   with slightly different initial conditions. This parameter expresses the
%   sensitivity to initial conditions. A positive Lyapunov
%   exponent is a signature of chaos.
%
%   Usage:
%   LYAPUNOV - calculate 50 days with a disturbance of 1E-15.
%   LYAPUNOV N NHORIZ - calculate N days with a disturbance of 1E-15. Use the 
%   first NHORIZ days to calulate the Lyapunov coefficient (time horizon).
%   LYAPUNOV N NHORIZ DISTURB = Disturb the initial conditions with DISTURB 
%   and run again.
%
%   See also lyapspect, lorenzmap, takens

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:27 $
function [res,nhoriz]=lyapunov(ndays, nhorizon, disturb)
global t g_Y g_t g_grind;
defaultdisturb=1E-8;
i_parcheck;
if nargin == 0
   ndays = g_grind.ndays;
else
   if ischar(ndays) && strncmpi(ndays,'-sep',4)
%      separate=1;
      ndays=g_grind.ndays;
   else 
      ndays = i_checkstr(ndays);
   end;
end;
if nargin <= 1
   nhorizon = ndays;
else
   nhorizon = i_checkstr(nhorizon);
end;
if nargin <= 2
 %  if isfield(g_grind.solver.opt,'AbsTol')
 %     disturb = g_grind.solver.opt.AbsTol*0.1;
 %  else
      disturb = defaultdisturb;
 %  end;
else
   disturb = i_checkstr(disturb);
end;
N0 = i_initvar;
oldtstep=g_grind.tstep;
g_grind.tstep=ndays;
try   
   i_ru(g_grind.odefile, t, ndays, N0, 0);
   Y = g_Y;
 % T = g_t;
   N0 = N0 + disturb;
   i_ru(g_grind.odefile, t, ndays, N0, 1);
   g_grind.tstep=oldtstep;
catch err
%   err=lasterror; 
   g_grind.tstep=oldtstep;
   rethrow(err);
end;
   
[a,b,dist,nhorizon]=calclyap(g_t,g_Y,Y,nhorizon);
H = i_makefig('lyap1');
plot(g_t, dist, g_t(1:nhorizon), exp(a * g_t(1:nhorizon) + b));
set(H, 'Name', 'Lyapunov plot');
if length(g_t)-nhorizon<10
   ch='> ';
else
   ch='= ';
end;
title(['\lambda max = ' num2str(a) ';  time horizon \tau ' ch num2str(nhorizon) ]);
xlabel('t');
ylabel('difference between 2 runs (log scale)');
set(gca, 'YScale', 'Log');
H = figure(i_figno('lyap2'));
plot(g_t, Y(:, 1), g_t, g_Y(:, 1));
title('Time plot with 2 slightly different initial settings');
set(H, 'Name', 'Lyapunov time plot');
xlabel('t');
ylabel(i_disptext(i_statevars_names(1)));
if nargout>0
   nhoriz=nhorizon;
   res=a;
end;

%%%% function calclyap
function [a,b,dist,nhorizon]=calclyap(gg_t,Y1,Y2,nhorizon)
dist = sqrt(sum((Y1 - Y2).^2,2));
%dists=abs(Y1-Y2);
%dists(dists<1E-16)=1E-16;
maxdist=max(dist);
if maxdist<1E-3*max(max(Y1));
   maxdist=1E-3*max(max(Y1));
end;
nhor=-1;
for tt = 1:length(gg_t);
   if abs(dist(tt)) < 1E-16
      warning('GRIND:lyapunov:nodiff','Small difference between runs: log of (almost) zero');
      dist(tt) = 1E-16;
   end;
   if (nhor==-1)&&(dist(tt)>maxdist/4)
      nhor=tt;
   end;
end;
if nhor~=-1
   nhorizon=nhor;
end;
n=nhorizon;
%tim = transpose([1:ndays]);
logd = log(dist(1:nhorizon));
tt1=gg_t(1:nhorizon);
sumx = sum(tt1);
sumy = sum(logd);
sumx2 = sum(tt1 .* tt1);
sumxy = sum(tt1 .*  logd);
a = (sumxy - sumx * sumy / n) / (sumx2 - sumx * sumx / n);
b = sumy / n - a * sumx / n;

