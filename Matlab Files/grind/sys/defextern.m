%DEFEXTERN   Define external variables. 
%   External variables are parameters that are variable in time. With the command setdata or 
%   loaddata the data are entered. This command can be used in the model definition, but can 
%   also be entered as a GRIND command to change a parameter in an external variable.
%
%   Usage:
%   DEFEXTERN VAR DEFAULT - Defines the variable VAR. DEFAULT is the default value used for 
%   simulation outside the scope of the data.
%   DEFEXTERN VAR DEFAULT MATRIX - You can also enter the data matrix directly.
%   DEFEXTERN VAR DEFAULT -cycle - option cycle reuses the data outside the range
%   DEFEXTERN VAR DEFAULT -nocycle - outside the range of the data the default value is used (default behaviour).
%   DEFEXTERN VAR DEFAULT -floor - do not interpolate within time steps.
%   DEFEXTERN -d - deactivate: external variable are (temporarily) considered to be parameters
%             so all data is neglected.
%   DEFEXTERN -a - reactivate external variables.       
%
%
%   See also model, definepars, defpermanent, setdata, loaddata

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:26 $
function g_result = defextern(name, default, data, opt1, opt2)
global g_grind g_Y;
cycle = 0;
tofloor = 0;
if nargin < 1
   if ~isfield(g_grind, 'externvars') || (isempty(g_grind.externvars))
      disp('No external variables defined');
   else
      disp('External variables (default values):');
      maxparlen=par('-maxparlen');
      for i = 1:length(g_grind.externvars)
         if isempty(g_grind.externvars{i}.data)
            fprintf(['%-' num2str(maxparlen) 's = %s  [No data] %s\n'], ...
               g_grind.externvars{i}.name, g_grind.externvars{i}.default,getoption(i));
         else
            fprintf(['%-' num2str(maxparlen) 's = %s  [datasize:%dx%d] %s\n'], ...
               g_grind.externvars{i}.name, g_grind.externvars{i}.default, ...
               size(g_grind.externvars{i}.data), getoption(i));
         end;
      end;
   end;
   return;
   % error('Cannot define external variable');
else
   if nargin == 1
      if strncmpi(name,'-d',2)
         %deactivate
         for i=1:length(g_grind.externvars)
            g_grind.externvars{i}.options.active=0;
            g_grind.pars=[g_grind.pars {g_grind.externvars{i}.name}];
         end;
         g_grind.lastsettings={};
         disp('All external variables are now parameters');
         return;
      elseif strncmpi(name,'-a',2)
         %activate
         for i=1:length(g_grind.externvars)
            g_grind.externvars{i}.options.active=1;
            j=find(strcmp(g_grind.externvars{i}.name,g_grind.pars));
            if ~isempty(j)
               g_grind.pars={g_grind.pars{1:j-1} g_grind.pars{j+1:end}};
            end;
         end;
         g_grind.lastsettings={};
         disp('All external variables are activated');
         return;
      end;
 %     if isempty(strfind(' ',name))&evalin('base',sprintf('exist(''%s'',''var'')',name))
 %        default=num2str(evalin('base',name));
 %     else
      default = '0';
 %     end
   end;
   f = strfind(name, ' ');
   if ~isempty(f)
      f2 = strfind(name, ';');
      if ~isempty(f2)
         name = name(1:f2(1) - 1);
      end;
      f2 = strfind(name, '%');
      if ~isempty(f2)
         name = name(1:f2(1) - 1);
      end;
      n = name;
      if (length(f) == 1)
         name = n(f(1) + 1:length(n));
         default = '0';
         data = NaN;
      else
         name = n(f(1) + 1:f(2) - 1);
         f1=strfind(name,'''');
         if length(f1) == 2
            name = name(f1(1) + 1:f1(2) - 1);
         end;
         if length(f) == 2
            default = n(f(2) + 1:length(n));
            f1=strfind(default,'''');
            if length(f1) == 2
               default = default(f1(1) + 1:f1(2) - 1);
            end;
            data = NaN;
         else
            default = n(f(2) + 1:f(3) - 1);
            if length(f) == 3
               opt1 = n(f(3) + 1:length(n));
               opt2 = ' ';
            else
               opt1 = n(f(3) + 1:f(4) - 1);
               opt2 = n(f(4) + 1:length(n));
            end;
            if strncmpi('-c',opt1, 2)||strncmpi('''-c',opt1, 2)||strncmpi('-c',opt2, 2)||strncmpi('''-c',opt2, 2)
               cycle = 1;
               data = NaN;
            elseif strncmpi('-n',opt1, 2)||strncmpi('''-n',opt1, 2)||strncmpi('-n',opt2, 2)||strncmpi('''-n',opt2, 2)
               cycle = 0;
               data = NaN;
            end;
            if strncmpi('-f',opt1, 2)||strncmpi('''-f',opt1, 2)||strncmpi('-f',opt2, 2)||strncmpi('''-f',opt2, 2)
               tofloor = 1;
               data = NaN;
            end;
         end;
      end;
   else
 %     if nargin < 2
 %       default = '0';
  %    end
      if iscell(default)
         default = default{1};
      end;
      if ~ischar(default)
         default = num2str(default);
      end;
      if (nargin == 3)
         if ischar(data)
            if strncmpi('-c',data, 2)
               cycle = 1;
               data = NaN;
            elseif strncmpi('-n',data,  2)
               cycle = 0;
               data = NaN;
            elseif strncmpi('-f',data,  2)
               tofloor = 1;
               data = NaN;
            else
               data = eval(data);
            end;
            
         end;
      else
         data = NaN;
      end;
      if (nargin == 4)
         if ischar(opt1)
            if strncmpi(opt1, '-c', 2)
               cycle = 1;
            elseif strncmpi(opt1, '-n', 2)
               cycle = 0;
            elseif strncmpi(opt1, '-f', 2)
               tofloor = 1;
            end;
         end;
      end;
      if (nargin == 5)
         if ischar(opt2)
            if strncmpi(opt2, '-c', 2)
               cycle = 1;
            elseif strncmpi(opt2, '-n', 2)
               cycle = 0;
            elseif strncmpi(opt2, '-f', 2)
               tofloor = 1;
            end;
         end;
      end;
      
   end;
end;
ivar = 0;
for i = 1:length(g_grind.externvars)
   if strcmp(name, g_grind.externvars{i}.name)
      ivar = i;
   end;
end;
if ivar == 0
   ivar = length(g_grind.externvars) + 1;
end;
ndx= ~strcmp(name,g_grind.pars);
ppars=g_grind.pars(ndx);
if length(g_grind.pars) > length(ppars)
   eval(['global ' name]);
   eval(['clear ' name]);
end;
g_grind.pars = ppars;
g_grind.externvars{ivar}.name = name;
g_grind.externvars{ivar}.default = default;
g_grind.externvars{ivar}.options.cycle = cycle;
g_grind.externvars{ivar}.options.tofloor = tofloor;
g_grind.externvars{ivar}.options.active =1;
if ~isnan(data)
   if ~isfield(g_grind.externvars{ivar}, 'data') || xor(isempty(g_grind.externvars{ivar}.data), isempty(data)) || (min(size(g_grind.externvars{ivar}.data) ~= size(data)) == 0) ...
         ||       (min(min(g_grind.externvars{ivar}.data == data)) == 0);
      g_Y = [];
   end;
   g_grind.externvars{ivar}.data = data;
else
   if ~isfield(g_grind.externvars{ivar}, 'data')
      g_grind.externvars{ivar}.data = [];
   end;
end;

modelchange = 1;
d1 = ['defextern ' name ];
d2=['defextern(''' name ];
for i = 1:length(g_grind.model)
   if strncmpi(d1, g_grind.model{i}, length(d1))
      d3 = [d1 ' ' default];
      if ~strncmpi(d3, g_grind.model{i}, length(d3));
         g_grind.model{i} = d3;
         i_mmodel;
         g_Y = [];
         if nargout == 1
            g_result= sprintf('%s=externvar(%d,%s,t);',name,ivar,default);
        %     [name '=externvar(' num2str(ivar) ', ' default ',t,g_grind.externvars{' num2str(ivar) '}.options);'];
         end;
         return;
      else
         modelchange = 0;
      end;
   end;
   if strncmpi(d2, g_grind.model{i}, length(d2))
      d3=[d2 ', ''' default ''''];
      if ~strncmpi(d3, g_grind.model{i}, length(d3));
         g_grind.model{i} = [d3 ');'];
         i_mmodel;
         g_Y = [];
         if nargout == 1
            g_result= sprintf('%s=externvar(%d,%s,t);',name,ivar,default);
            % g_result=[name '=externvar(g_grind.externvars{' num2str(ivar) '}.data, ' default ',t,g_grind.externvars{' num2str(ivar) '}.options);'];
         end;
         return;
      else
         modelchange = 0;
      end;
   end;
end;
if modelchange
   g_grind.model = [{[d1 ' ' default]}, g_grind.model];
   starting=isfield(g_grind,'starting');
   i_mmodel;
   if starting
       g_grind.starting=1;
   end;
   g_Y = [];
end;
if nargout == 1
%   g_result=[name '=externvar(g_grind.externvars{' num2str(ivar) '}.data, ' default ',t,g_grind.externvars{' num2str(ivar) '}.options);'];
  g_result= sprintf('%s=externvar(%d,%s,t);',name,ivar,default);
end;
function s=getoption(ivar)
global g_grind;
s='';
if g_grind.externvars{ivar}.options.cycle 
   s=[s '-cycle '];
end;
if g_grind.externvars{ivar}.options.tofloor 
   s=[s '-floor '];
end;
