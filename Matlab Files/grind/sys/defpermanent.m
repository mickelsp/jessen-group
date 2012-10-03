%DEFPERMANENT   Define permanent variables. 
%   "Permanent variables" are variables that may change during a run. They are not state variables
%   as they are not part of the differential or difference equation. They offer for instance an easy
%   way to monitor the maximum of a state variable. The initial value can be set in the same 
%   way as with state variables. You can plot the permanent variables in time plots etc. Their 
%   values are stored in the g_permanent global variable.
%
%   Usage:
%   DEFPERMANENT VAR - Defines the variable VAR. DEFAULT is the default value used for 
%   simulation outside the scope of the data.
%   DEFPERMANENT VAR VALUE - You can also enter the data matrix directly.
%
%
%   See also model, definepars, defextern, outfun, time
%

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:26 $
function res = defpermanent(avar, avalue)
global g_permanent g_grind g_t;
if nargin  > 0
   if strcmpi(avar, '-initiate')
      %initiate before run
      
      %   in  curr_ode.m:
      %   remove defpermanent lines
      %   replace "avar" by "g_grind.permanent{avarno}.currvalue" with i_changevar
      %   update all at the end
      %   i_updatepermanent(at);
      %
      tspan = avalue;
      offs = 1;
      for i = 1:length(g_grind.permanent)
         g_grind.permanent{i}.currvalue = evalin('base', g_grind.permanent{i}.name);
         if isempty(g_grind.permanent{i}.currvalue)
            g_grind.permanent{i}.dims = [1,1];
         else
            g_grind.permanent{i}.dims = size(g_grind.permanent{i}.currvalue);
         end;
         g_grind.permanent{i}.from = offs;
         offs = offs + prod(g_grind.permanent{i}.dims);
         g_grind.permanent{i}.to = offs - 1;
      end;
      if length(tspan) == 2
        tspan(2)=tspan(1)+tspan(2);
        if isnan(g_grind.tstep)
            tspan = tspan(1):tspan(2);
         else
            tspan = tspan(1):(tspan(2) - tspan(1)) / g_grind.tstep:tspan(2);
         end;
      end;
      g_permanent.t = tspan(:);
      g_permanent.Y = zeros(length(g_permanent.t), offs - 1) + NaN;
      g_permanent.nextt = 1;
      g_permanent.active = 1;
      g_permanent.lasti = 1;
      g_permanent.lastt(g_permanent.lasti) = tspan(1) - 1;
      g_permanent.lastY=cell(10,1);
      i_updatepermanent(tspan(1)); %set t = 0 to initial values
   elseif strcmpi(avar, '-updatevars')
      for i = 1:length(g_grind.permanent)
         g_grind.permanent{i}.currvalue = evalin('base', g_grind.permanent{i}.name);
      end;
   elseif strcmp(avar, '-deactivate')
      if g_permanent.active&&(g_permanent.nextt <= length(g_permanent.t))
         i_updatepermanent;
         g_permanent.nextt = g_permanent.nextt - 1;
      end;
      g_permanent.active = 0;
   elseif strcmp(avar, '-activate')
      g_permanent.active = 1;
      if ~isempty(avalue)
         defpermanent('-s', avalue);
      end;
   elseif strncmpi(avar, '-g', 2)
      %output
      defpermanent('-deactivate');
      if isempty(g_permanent) || isempty(g_permanent.Y) || (~isempty(avalue)&&(avalue > length(g_grind.permanent)))
         res = [];
      else
         if isempty(avalue)
            res = g_permanent.Y;
         else   
            res = g_permanent.Y(:, g_grind.permanent{avalue}.from:g_grind.permanent{avalue}.to);
         end;
         if length(g_t) ~= length(g_permanent.t)
            res = interp1(g_permanent.t, res, g_t);
         end;
      end
   elseif strncmpi(avar, '-l', 2) %-list
      disp('');
      maxparlen = par('-maxparlen');
      if ~isfield(g_grind,'permanent')||isempty(g_grind.permanent)
         disp('No permanent variables');
      elseif isempty(g_permanent.Y)
         disp('Initial values of permanent variables:')
         s=['%-' num2str(maxparlen) 's = %0.6g\n'];
         for i = 1:length(g_grind.permanent)
            p = evalin('base', g_grind.permanent{i}.name);
            fprintf(s, g_grind.permanent{i}.name, p(:));
         end;
      else
         disp('Initial and final values of permanent variables in last run:')
         s=['%-' num2str(maxparlen) 's = %0.5g / %0.5g\n'];
         for i = 1:length(g_grind.permanent)
            p = evalin('base', g_grind.permanent{i}.name);
            fprintf(s, g_grind.permanent{i}.name, p(:), defpermanent('-e', i));
         end;
      end;
   elseif strncmpi(avar, '-p', 2)
      %get the values at t=at
      defpermanent('-deactivate');
      if isempty(g_grind.permanent)
         res = [];
      elseif (nargin==1)||isempty(avalue) %defpermanent('-p') gets the initial values
         if isfield(g_grind.permanent{end},'to')
             res = zeros(g_grind.permanent{end}.to, 1);
             for i = 1:length(g_grind.permanent)
                res(g_grind.permanent{i}.from:g_grind.permanent{i}.to) = g_grind.permanent{i}.currvalue(:);
              end;
         else
             res=zeros(size(g_grind.permanent));
             for i = 1:length(g_grind.permanent)
                res(i) = g_grind.permanent{i}.currvalue(:);
              end;
         end;
      elseif ~(isempty(g_permanent) || isempty(g_permanent.Y))
         at = avalue;
         f=find(g_permanent.t >= at);
         if isempty(f)
            at = length(g_permanent);
         else
            at = f(1);
            if at == g_permanent.nextt %add the last values without increasing nextt
               for i = 1:length(g_grind.permanent)
                  g_permanent.Y(at, g_grind.permanent{i}.from:g_grind.permanent{i}.to) = g_grind.permanent{i}.currvalue(:)';
               end;
            end;
         end;
         if at >= 1
            res = g_permanent.Y(at, :)';
         else
            res = [];
         end;
      else
           res=[];
      end
   elseif strncmpi(avar, '-s', 2)
      %set the current values to avalue
      if isempty(avalue)
         avalue = defpermanent('-p', length(g_permanent.t));
      end;
      if  ~isempty(avalue)
         for i = 1:length(g_grind.permanent)
            p = avalue(g_grind.permanent{i}.from:g_grind.permanent{i}.to);
            p = reshape(p, g_grind.permanent{i}.dims);
            assignin('base', g_grind.permanent{i}.name, p);
         end;
      end;
   elseif strncmpi(avar, '-e', 2)
      %last value
      defpermanent('-deactivate');
      if isempty(g_permanent) || isempty(g_permanent.Y)|| avalue > length(g_grind.permanent)
         res = [];
      else
         res = g_permanent.Y(end, g_grind.permanent{avalue}.from:g_grind.permanent{avalue}.to);
      end
   elseif strcmp(avar, '-i_mmodel')
      %define permanent variable (internally used by i_mmodel)
      %defpermanent(avar)
      avar = strtrim(avalue);
      f =  strfind(avar, ';');
      if ~isempty(f)
         avar = avar(1:f(1) - 1);
      end;
      f =  strfind(avar, ' ');
      if ~isempty(f)
         avalue = str2num(avar(f(1) + 1:end)); %#ok
         avalue =  i_checkstr(avalue);
         avar = avar(1:f(1) - 1);
         evalin('base',sprintf('global %s g_permanent',avar));
         assignin('base', avar, avalue);
      else
         evalin('base',sprintf('global %s g_permanent',avar));
         assignin('base', avar, NaN);
      end;
      if ~isempty(g_grind.permanent)
         no = length(g_grind.permanent) + 1;
      else
         no = 1;
      end;
      g_permanent.Y = [];
      g_permanent.t = [];
      g_permanent.nextt = 1;
      g_permanent.active = 0;
      g_permanent.lastt = -9999+zeros(10,1);
      g_permanent.lastYs = cell(10,1);
      g_permanent.lasti = 1;
      g_grind.permanent{no}.name = avar;
      if isempty(avalue)
         g_grind.permanent{no}.currvalue = 0;
      else
         g_grind.permanent{no}.currvalue = avalue;
      end
      g_grind.permanent{no}.initiate=1;
%       ppars = {};
%        for i = 1:length(g_grind.pars)
%          if ~strcmp(g_grind.pars{i}, avar)
%            ppars = [ppars g_grind.pars(i)];
%          end;
%         end;
      g_grind.pars = g_grind.pars(~strcmp(avar,g_grind.pars));
   else
      i_parcheck;
      for i = 1:length(g_grind.permanent)
         if strcmp(avar, g_grind.permanent{i}.name)
            fprintf('"%s" is already defined as permanent variable\n', avar);
            return;
         end
      end;
      ButtonName = questdlg('Defining a permanent variable requirs a reset of the current run, OK to continue?', ...
         sprintf('defpermanent %s',avar),'OK','Cancel','OK') ;
      if strcmp(ButtonName, 'OK')
         disp('OK');
         if (nargin==1) && evalin('base',sprintf('exist(''%s'',''var'')',avar))
            avalue = evalin('base', avar);
         else
            avalue = [];
         end;
         finishgrind;
         g_grind.model=[{strtrim(sprintf('defpermanent %s %g',avar,avalue))}, g_grind.model];
         i_mmodel(g_grind.model, g_grind.commands, g_grind.inifile, g_grind.scheme);
      end;
   end;
end;
