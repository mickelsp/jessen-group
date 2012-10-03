%EULER   Euler integration
%   Solve an odefunct using simple Euler integration.
%   Use this integration method only with very small
%   time steps or to study numerical problems with differential equations.
%
%   See also solver, rk4, ode45, ode23, ode113, ode15s, ode23s, 
%   ode23t, ode23tb, odeset, odeget

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:26 $
function [tout, yout] = euler(odefunct, tspan, y0, options)
if (nargin < 4) || isempty(options.MaxStep)
   %if there are no valid options use default stepsize
   delta = 0.1;
else
   %use the option MaxStep as step size
   delta = options.MaxStep;
end

% Test that tspan is internally consistent.
tspan = tspan(:);
ntspan = length(tspan);
if ntspan == 1
   t0 = 0;
   next = 1;
else
   t0 = tspan(1);
   next = 2;
end
tfinal = tspan(ntspan);
if t0 == tfinal
   error('GRIND:euler:tpan','The last entry in tspan must be different from the first entry.');
end
tdir = sign(tfinal - t0);
if any(tdir * (tspan(2:ntspan) - tspan(1:ntspan-1)) <= 0)
   error('GRIND:euler:tspan','The entries in tspan must strictly increase or decrease.');
end
t = t0;
y = y0(:);
neq = length(y);

% Set the output flag.

outflag = ntspan > 2;                          % output only at tspan points

% Allocate memory if we're generating output.
delta = delta * tdir;

if nargout > 0
   if ntspan > 2                         % output only at tspan points
      tout = zeros(ntspan, 1);
      yout = zeros(ntspan, neq);
%      nsteps=round(tdir*(tfinal-t0)/delta)+1;
   else
      tout = transpose(t0:delta:tfinal);
%      nsteps=size(tout,1);
      yout = zeros(size(tout, 1),neq);
   end
   nout = 1;
   tout(nout) = t;
   yout(nout, :) = y.';
end
%
%MAIN LOOP
%evaluate the odefunction for the next time steps
%fold=[];
running=1;
while running;
   f=feval(odefunct, t, y);
%  simple Euler
   ynew = y + f .*  delta;
%  improved Euler (trapesium)
%  if isempty(fold)
%    ynew = y + f .*  delta;
%  else
%    ynew = y + 0.5*(fold+f) .*  delta;
%  end;
%  fold=f;
   tnew = t + delta;
   if tnew+0.01*delta>=tfinal
      running=0;
      tnew=tfinal;
   end;
   if ~outflag                 % computed points, no refinement
      nout = nout + 1;
      tout(nout) = tnew;
      yout(nout, :) = ynew.';
   elseif (tdir * (tnew - tspan(next)) >= 0) % at tspan, tspan assumed to be larger than delta
       nout = nout + 1;
       tout(nout) = tnew;
       yout(nout, :) = ynew.';
       next = next + 1;
   end;
   y = ynew;
   t = tnew;
end;   
if nout<length(tout)
   tout=tout(1:nout);
   yout=yout(1:nout,:);
end;

   
