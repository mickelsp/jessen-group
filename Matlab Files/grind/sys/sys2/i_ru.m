function i_ru(odefile, at, ndays, N0,updatesettings)
global g_t g_Y g_grind g_func t g_data;
if ~isempty(g_grind.onstart.funs)
   onstart;
end;
if ~isempty(g_grind.permanent)
   defpermanent('-initiate',[at,ndays]);
end;
if updatesettings
   g_func = [];
   g_grind.lastsettings = i_getsettings(N0);
   implicitdisperse;
end;
if g_grind.solver.reset_at_data
   if isfield(g_data, 'obs')
      tt = g_data.t;
      if at < tt(1)
         tt = [at; tt];
      end;
      g_tt = at;
      g_YY = N0';
      oldtstep = g_grind.tstep;
      try
         g_grind.solver.reset_at_data = 0;
         g_grind.tstep = 3;
         for i = 1:length(tt) - 1
            N0_obs = g_data.obs(i, :)'; %reset to the observation
            N0_obs(isnan(N0_obs)) = N0(isnan(N0_obs));
            i_ru(g_grind.odefile, tt(i),tt(i + 1) - tt(i), N0_obs,0);
            g_tt = [g_tt; g_t(1, :) + 0.01; g_t(end, :)];
            g_YY = [g_YY; g_Y(1, :); g_Y(end, :)];
            N0 = g_Y(end, :)';
         end;
         if ndays - at > tt(end)
            g_grind.tstep = oldtstep - length(g_tt);
            i_ru(g_grind.odefile, tt(end),ndays - tt(end), N0_obs,0);
            g_tt = [g_tt; g_t];
            g_YY = [g_YY; g_Y];
         end;
         g_Y = g_YY;
         g_t = g_tt;
         g_grind.tstep = oldtstep;
         g_grind.solver.reset_at_data = 1;
      catch %#ok
         g_grind.solver.reset_at_data = 1;
         g_grind.tstep = oldtstep;
      end;
      return;
   end;
end
if g_grind.solver.addmode
   oldt = g_t;
   oldy = g_Y;
   if ~isempty(g_t)
      if g_grind.solver.backwards
         at = -g_t(1);
         N0 = g_Y(1, :)';
      else
         at = g_t(length(g_t)) + 0.001;
         N0 = g_Y(size(g_Y, 1), :)';
      end;
   end;
end;
if g_grind.solver.backwards
   if g_grind.solver.isdiffer
      odefile = 'i_backdiffer';
   else
      if g_grind.solver.haslag
         g_grind.solver.name = 'ode23';
      end;
      odefile = 'i_back';
   end;
end;
try
   if isfield(g_grind, 'boxcar')
      for i = 1:length(g_grind.boxcar.trains)
         g_grind.boxcar.trains{i}.gcycl = 0;
      end;
   end;
   if g_grind.solver.hasevents
      settings = i_getsettings(N0);
      N00=N0;
      [N0,NP0] = i_initvar;

     try
         nextevent=min(ndays,i_setevent('runnext','',-1));
         t0 = at;
         te = at + ndays;
         g_t = [];
         g_Y = [];
         oldt=t;
         while t0 < te
            if t0 <= nextevent
               [N0,NP] = i_initvar;
               if ~isnan(g_grind.tstep)
                  %create equally spaced output
                  ts = t0:ndays /g_grind.tstep:nextevent;
                  if ts(end) < nextevent
                     ts = [ts, nextevent];
                  end;
                  if length(ts) < 5 %at least 5 steps
                     ts = t0:(nextevent - t0) / 5:nextevent;
                  end;
               else
                  ts = [t0; nextevent];
               end;
               if (length(ts) > 1) && (ts(1) < ts(end))
                  if ~isempty(g_grind.permanent) 
                     defpermanent('-activate',NP);
                  end;
                  i_keep(N00) %needed if you want to access the initial value from the program (val)
                  [tt, YY] = feval(str2func(g_grind.solver.name), str2func(odefile), ts, N0, g_grind.solver.opt);
                  g_t = [g_t; tt];
                  g_Y = [g_Y; YY];
                  ke; %now the event can change the values of the state variables
                  t=t0; %#ok
               end;
               t0 = nextevent + 1E-10;
            end;
            nextevent=i_setevent('runnext','',nextevent);
            if ~isempty(g_grind.permanent)
               defpermanent('-updatevars');
            end; 
            if isnan(nextevent) || (nextevent > te)
               nextevent = te;
            end;
         end;
         if ~isnan(g_grind.tstep)
            ts = (at:ndays / g_grind.tstep:te)';
            g_Y(end,:) = i_initvar';
            g_Y = interp1(g_t, g_Y, ts);
            g_t = ts;
         end;
         i_setsettings(settings);
         if ~isempty(g_grind.permanent)
           defpermanent('-s',NP0);
         end;
         t=oldt;
     catch err
%          err=lasterror;
          i_setsettings(settings);
          if ~isempty(g_grind.permanent)
            defpermanent('-s',NP0);
          end;
          t=oldt; %#ok
         rethrow(err);
      end;
   else
      if ~isnan(g_grind.tstep)
        ts = (at:ndays / g_grind.tstep:at + ndays)';
      else
         ts = [at; at + ndays];
      end;
%       if g_grind.version.isoctave
%          g_grind.solver.opt.MaxStep = 0.2; %necessary??
%          g_grind.solver.opt.OutputFcn = [];
%       end;
      %#function ode45, euler, i_differ
      [g_t, g_Y] = feval(str2func(g_grind.solver.name), str2func(odefile), ts, N0, g_grind.solver.opt);
   end;
catch err
%    err =lasterror;
     g_grind.lastsettings.initvar = [];
    rethrow(err);
end;
if g_grind.solver.backwards
   if g_grind.solver.haslag
      g_grind.solver.name = 'ddesol';
   end;
   g_t = -flipud(g_t);
   g_Y = flipud(g_Y);
end;
if g_grind.solver.addmode
   if g_grind.solver.backwards
      g_t = [g_t; oldt];
      g_Y = [g_Y; oldy];
   else
      g_t = [oldt; g_t];
      g_Y = [oldy; g_Y];
   end;
end;


