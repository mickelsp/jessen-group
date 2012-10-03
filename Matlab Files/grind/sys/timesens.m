%TIMESENS   Performs numerical sensitivity analysis of a parameter in time
%  This algorithm does two precise runs with a slightly different parameter
%  or initial condition. It computes the derivative during the run. The plots are
%  added to a time plot if possible. Optionally the elasticities (dV/dp*p/V) are
%  plotted.
%
%  Usage:
%  TIMESENS PAR - adds PAR in the time plots
%  TIMESENS PAR PAR2 PAR3 etc - analyses PAR1-3
%  TIMESENS -ALL - adds all parameters and state variables
%  TIMESENS -C - clear all data
%  TIMESENS -U - update all
%  TIMESENS -D 1E-10 sets disturbance to 1E-10 (default=1E-8)
%  TIMESENS -D PAR 1E-10 - sets only disturbance for PAR.
%  TIMESENS -E .... - plots elasticities instead of derivatives.
%
%  See also time, simtime, out

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:29 $
function [result] = timesens(varargin)
global g_Y g_grind g_sens;
evalin('base','global g_sens');
defaultdisturb = 1E-8;
i_parcheck;
ndx = strcmpi('-e', varargin);
if any(ndx)
   g_sens.elasticity = 1;
   varargin = varargin(~ndx);
else
   g_sens.elasticity = 0;
end;
if isempty(varargin)
   answer=inputdlg({'Enter parameter:'},'timesens',1);
   if isempty(answer) && ~isempty(answer{1})
      return;
   end;
   varargin{1} = answer{1};
end;
par = varargin{1};
if ischar(par)
   if par(1) == '-'
      if strcmpi(par, '-u')
         doupdate;
         return;
      elseif strcmpi(par, '-all')
         varargin = [g_grind.pars g_grind.statevars.names];
      elseif strcmpi(par, '-d')
         if length(varargin) == 2
            disturb = i_checkstr(varargin{2});
            for i = 1:length(g_sens.pars)
               g_sens.disturb(i) = disturb;
            end;
         else
            i = findadd(varargin{2});
            if i <= length(g_sens.pars)
               g_sens.disturb(i) = i_checkstr(varargin{3});
            end;
         end
         return;
      elseif strcmpi(par, '-c')
         if ~isempty(g_sens)
            for i = 1:length(g_grind.timevars)
               for k = 1:length(g_sens.pars)
                  for j = 1:length(g_grind.statevars.names)
                     out('-silent',sprintf('-%d',i),'-remove',sprintf('timesens(''d(%s)/d(%s)'')',...
                        g_grind.statevars.names{j}, g_sens.pars{k} ));
                  end;
               end;
            end;
         end;
         g_sens = [];
         out('-silent','-cleantimevars')
         disp('Cleared all');
         return;
      end;
   end;
   f2 = strfind(par, ')/d(');
   if isempty(f2)
      f2 = strfind(par, ')./d(');
   end;
   if ~isempty(f2)   % d(par) / d(var) or  d(par(i)) / d(var(i))
      f1 = strfind(par, 'd(');
      f3 = strfind(par, ')');
      var = par(f1(1) + 2:f2(1) - 1);
      par1 = par(f1(2) + 2:f3(end) - 1);
      ivar = i_getno(var);
      if ~ivar.isvar
         error('GRIND:timesens:NoStatevars','%s is not a state variable',var);
      end;
      ipar = findadd(par1);
      if ~isfield(g_sens, 'Ys') || length(g_sens.Ys) < ipar
         addpar(par1, defaultdisturb);
         g_sens.lastsettings = [];
      end;
      if settingschanged
         doupdate;
      end;
      if ~isfield(g_sens, 'disturb') || (length(g_sens.disturb) < ipar) || (g_sens.disturb(ipar)==0)
         g_sens.disturb(ipar) = defaultdisturb;
      end;
      yy = g_sens.Ys{ipar};
      result = (yy(:,ivar.no) - g_Y(:, ivar.no)) / g_sens.disturb(ipar);
      return;
   end;
end;
npar = 0;
if g_sens.elasticity
   s = 'Added elasticities of';
else
   s = 'Added sensitivities of';
end;
if isfield(g_sens,'pars')
    oldn=length(g_sens.pars);
else
    oldn=0;
end
for j = 1:length(varargin)
   if addpar(varargin{j}, defaultdisturb)
      npar = npar + 1;
   end;
end;
if npar > 0
   if oldn==length(g_sens.pars)
      disp('Parameter(s) were already selected');
      return;
   end;
   s = sprintf('%s ', s, g_sens.pars{oldn+1:end});
   disp(s);
   disp('Run <a href="matlab: time">time</a> to see the results');
   if strcmp(g_grind.solver.name, 'ode45')
      if g_grind.solver.opt.AbsTol > 1e-10 || g_grind.solver.opt.RelTol > 1e-10
         disp('timesens has set solver tolerances to values 1E-11');
         solver('ode45', 1E-11, 1E-11);
      end;
   end;
end;

function ok = addpar(par, disturb)
global g_grind g_sens;
no = i_getno(par);
ok = 1;
if isempty(no.no)
   fprintf('unknown parameter %s\n', par);
   ok = 0;
   return;
end;
siz= evalin('base',sprintf('size(%s)',par));
if prod(siz) > 1
   if (siz(2) == 1) || (siz(1)==1)
      for j = 1:siz(1)*siz(2)
         ok = addpar(sprintf('%s(%d)', par, j),disturb);
      end
   else
      for i = 1:siz(1)
         for j = 1:siz(2)
            ok = addpar(sprintf('%s(%d,%d)',par,i,j),disturb);
         end
      end;
   end;
else
   i = findadd(par);
   parexist = ~isempty(g_sens) && isfield(g_sens, 'pars');
   if parexist
      parexist=(i <= length(g_sens.pars));
   end;
   g_sens.pars{i} = par;
   g_sens.disturb(i) = disturb;
   if ~parexist
      g_sens.lastsettings = [];
      len = length(g_grind.timevars) + 1;
      for j = 1:g_grind.statevars.dim
         avar = i_statevars_names(j);
         if ~g_sens.elasticity
            out('-silent',sprintf('-%d',len),'-add',sprintf('timesens(''d(%s)/d(%s)'')',...
               avar, g_sens.pars{i}));
         else
            out('-silent',sprintf('-%d',len),'-add',sprintf('timesens(''d(%s)/d(%s)'')*%s/%s',...
               avar, g_sens.pars{i}, g_sens.pars{i},avar));
         end;
      end;
   end;
end;
function res = settingschanged
global g_sens g_t;
if isempty(g_sens) || ~isfield(g_sens, 'lastsettings') || isempty(g_sens.lastsettings)
   res = 1;
   return;
end;
settings = i_getsettings;
if ~isempty(g_t);
   settings(1) = g_t(end);
end;
res = isdifferent(settings, g_sens.lastsettings);

function res = isdifferent(A, B)
res=~min(size(A) == size(B));
if isempty(B) && isempty(A)
   res = 0;
elseif ~res
   compar = A == B;
   res = ~min(compar);
   if res
      res = ~min(compar + isnan(A) .* isnan(B));
   end;
end;

function doupdate
global g_Y g_t t g_grind g_sens;
disp('updating timesens...');
g_sens.lastsettings = i_getsettings;
g_sens.lastsettings(1) = g_t(end);
oldY = g_Y;
oldt = g_t;
try
   for l_i = 1:length(g_sens.pars)
      oldpar = evalin('base', g_sens.pars{l_i});
      multassignin('base',  g_sens.pars{l_i}, oldpar + g_sens.disturb(l_i));
      N0 = i_initvar;
      i_ru(g_grind.odefile, t, g_t(end), N0, 0);
      multassignin('base', g_sens.pars{l_i}, oldpar);
      g_sens.Ys{l_i} = interp1(g_t, g_Y, oldt);
   end;
   g_Y = oldY;
   g_t = oldt;
catch err
   %   err=lasterror;
   g_Y = oldY; %#ok
   g_t = oldt; %#ok
   rethrow(err);
end;

function i = findadd(par)
global g_sens;
i = 1;
if ~isempty(g_sens) && isfield(g_sens, 'pars')
   while i <= length(g_sens.pars) && ~strcmp(g_sens.pars{i}, par)
      i = i + 1;
   end;
end;

