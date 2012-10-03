%OUT   Select variables for time plots
%   Specify the variables/functions for the time plots. You can define different
%   time plots that are updated simultaneously.
%
%   Usage: 
%   OUT - opens a dialog box to add variables to time plots 
%   OUT X1 X2 X3 - plots the variables X1,X2 and X3 (replaces current
%   list). X1 may be a state variable, a parameter, an intermediate variable
%   (see FUNCS), or an MATLAB expression with these variables; 
%   OUT -OPT X4 X5 - uses option -OPT
%
%   Valid options are:
%     -1, -2, -VALUE, define the variable/functions for time figure number VALUE 
%      (default VALUE=1, can be in combination with other options, but this 
%      should then be the first option)
%     -add, adds the variables/functions X4 and X5 to the current list
%     -replace, replaces the current list with X4 and X5
%     -remove, removes X4 and X5 from the current list
%     -time, replaces the time axis with the function X4 (X5 is neglected)
%     -clear, clear all time plots. To remove one plot, use in combination with the first option.
%
%   Shorthands (can be used in combinations and with previous options)
%     -defaults, all state variables (including data) and external variables.
%     -data, all state variables for with data is defined.
%     -extern, all external variables.
%     -sum, sets the sum of all vector state variables (works also with all other options of 
%      outf, e.g. -mean -max -min -median).
%
%   Examples:
%   OUT A  - only the state variable A on the axis
%   OUT -add cons - add the intermediate variable cons
%   OUT -2 A - make a second plot with state variable A
%   OUT -3 -mean -sum -max - mean sum and maximum of vector state variables in plot 3.
%   OUT -add A/K - add an expression with parameters and state variables
%   OUT -time mod(t,365) - replaces the time axis with daynumber
%
%   See also time, funcs, outf, timesens
%

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:28 $
function out(varargin)
global g_grind;
if isempty(g_grind)
    error('GRIND:out:NoModel','No model selected');
end;
silent = 0;
No = 1;
hasNo = 0;
if (nargin > 0)
   if iscell(varargin{1})
      varargin = varargin{1};
   end;
   if strcmp(varargin{1}, '-silent')
      varargin = varargin(2:length(varargin));
      if isempty(varargin)
         i_outdlg;
      end;
      silent = 1;
   end;
   if ~isempty(strfind(varargin{1}, '{'))
      varargin = eval(char(varargin{1}));
   end;
   if strncmp(varargin{1}, '-',1) && isnumstring(varargin{1})
      No = str2double(varargin{1});
      if isempty(No)
         No = 1;
      else
         hasNo = 1;
         No = -No;
         varargin = varargin(2:end);
         if isempty(varargin)
            i_outdlg('init', No);
            return;
         end;
      end;
   end;
   if strncmpi(varargin(1), '-ch',3)
      if hasNo
         g_grind.timevars{No} = checklst(g_grind.timevars{No});
      else
         for iNo = length(g_grind.timevars):-1:1
            if isempty(g_grind.timevars{iNo})
               g_grind.timevars =  {g_grind.timevars{1:iNo - 1} g_grind.timevars{iNo + 1:end}};
               g_grind.outt =  {g_grind.outt{1:iNo - 1} g_grind.outt{iNo + 1:end}};
            end;
         end;
         for iNo = 1:length(g_grind.timevars)
            g_grind.timevars{iNo} = checklst(g_grind.timevars{iNo});
         end;
      end;
   elseif strcmpi(varargin(1), '-cleantimevars')
      cleantimevars;
   elseif strncmpi(varargin(1), '-cl',3)
      if hasNo
         g_grind.timevars{No} = {};
         g_grind.outt{No} = '';
         if ~silent
            fprintf('Time plot %d cleared\n', No);
         end;
      else
         g_grind.timevars = {};
         g_grind.outt = {};
      end;
   elseif strcmpi(varargin{1}, '-remove')
      if length(varargin) > 1 && strcmpi(varargin{2}, '-data')
         varargin{2} = '-removedata';
      end;
      ntime = size(g_grind.timevars{No}, 2);
      varlist = getvarlist(No, varargin{2:end});
      for j = 1:length(varlist)
         i = 1;
         while (i <= ntime) && ~strcmp(g_grind.timevars{No}{i}, varlist{j})
            i = i + 1;
         end;
         [g_grind.timevars{No},ntime] = removei(g_grind.timevars{No},i, ntime);
      end;
   elseif strcmp(varargin{1}, '-removeextern')
      out('-silent',num2str(-No),'-remove','-extern')
   elseif strcmp(varargin{1}, '-removedata')
      ntime = size(g_grind.timevars{No}, 2);
      for k = length(g_grind.timevars{No}):-1:1
         if length(g_grind.timevars{No}{k}) >= 9
            if strcmpi(g_grind.timevars{No}{k}(1:9), 'observed ')
               [g_grind.timevars{No},ntime] = removei(g_grind.timevars{No},k, ntime);
            end;
         end;
      end;
   elseif strcmp(varargin{1}, '1')||strcmpi(varargin{1}, '-add')
      addtovars(getvarlist(No,varargin(2:end)), No);
   elseif strcmp(varargin{1}, '0')||strcmpi(varargin{1}, '-replace')
      g_grind.timevars{No} = getvarlist(No, varargin(2:end));
   elseif strcmpi(varargin{1}, '-time')
      if length(varargin)==1 
         varargin={'-time','t'};
       end;
      if strncmp(varargin{2}, '-',1) && isnumstring(varargin{2})
         No = str2double(varargin{2});
         if isempty(No)
            No = 1;
         else
            No = -No;
            varargin = {'-time', varargin{3:end}};
            if length(varargin) == 1
               varargin={'-time','t'};
            end;
         end;
      end;
      g_grind.outt{No} = varargin{2:end};
   elseif ~strcmp(varargin{1}, '?')
      g_grind.timevars{No} = getvarlist(No, varargin);
   end;
else
   %  outdlg;
   i_outdlg;
   return;
end;
if g_grind.statevars.vector
   done = 0;
   timvars = {};
   if ~isempty(g_grind.timevars)
      for i = size(g_grind.timevars{No}, 2):-1:1
         if (i <= length(g_grind.timevars{No})) && isempty(strfind(g_grind.timevars{No}{i},'('))
            for j = 1:size(g_grind.statevars.vectnames, 2)
               if (i <= length(g_grind.timevars{No})) && strcmp(g_grind.timevars{No}{i}, g_grind.statevars.vectnames{j})
                  if g_grind.statevars.dims{j}.dim2 == 1
                     done = 1;
                     g_grind.timevars{No} =  removei(g_grind.timevars{No}, i, length(g_grind.timevars{No}));
                     for k = 1:g_grind.statevars.dims{j}.dim1
                        addtimvars = {sprintf('%s(%d)', g_grind.statevars.vectnames{j}, k)};
                        if isempty(timvars)
                           timvars = addtimvars;
                        else
                           timvars = [timvars, addtimvars];
                        end;
                     end;
                  end;
               end
            end
         end
      end;
   end
   if done
      g_grind.timevars{No} = [g_grind.timevars{No}, timvars];
   end;
end;
if (length(varargin) ~= 1) || ~strcmp(varargin{1}, '?')
   if No <= length(g_grind.timevars)
      for i = 1:length(g_grind.timevars{No})
         g_grind.timevars{No}{i} = outf('changeshortcut', g_grind.timevars{No}{i});
         if No > length(g_grind.outt) || isempty(g_grind.outt{No});
            g_grind.outt{No} = 't';
         end;
      end;
   end;
end;
if ~silent
   if isempty(g_grind.timevars)
      disp('No time variables');
   end;
   for iNo = 1:length(g_grind.timevars)
      if ~isempty(g_grind.timevars{iNo}) && (~hasNo || (iNo == No))
         if length(g_grind.timevars{iNo}) > 10
            s1 = i_cell2str({g_grind.timevars{iNo}{1:9} '...' g_grind.timevars{iNo}{end}});
         else
            s1 = i_cell2str(g_grind.timevars{iNo});
         end;
         s = sprintf('Time (%s) plot %d:  %s',g_grind.outt{iNo}, iNo, s1);
         disp(s);
      end;
   end;
end;

function reslist = checklst(varlist)
global g_data;
k = 1;
reslist = cell(1, length(varlist));
for i = 1:length(varlist)
   OK2add = ~isempty(strtrim(varlist{i}));
   if OK2add
      % check for 'observed ' without external data
      if strncmpi(varlist{i}, 'observed ',9)
         if ~isempty(g_data) && ~isempty(g_data.obs)
            avar = strtrim(varlist{i}(9:end));
            iX = i_getno(avar);
            OK2add = iX.isvar & ~min(isnan(g_data.obs(:, iX.no)));
         else
            OK2add = 0;
         end;
      end;
      %check for double entries
      if OK2add
         for j = 1:i - 1
            if strcmp(varlist{j}, varlist{i})
               OK2add = 0;
               break;
            end;
         end;
      end;
   end;
   if OK2add
      reslist{k} = varlist{i};
      k = k + 1;
   end;
end;
if k > 1
   reslist = reslist(1:k - 1);
else
   reslist = {};
end;

function  varlist1 = getvarlist(No, varlist)
global g_grind g_data;
varlist1 = cell(1, length(varlist));
k = 1;
if nargin > 0
   if ~iscell(varlist)
      varlist = {varlist};
   end;
   for i = 1:length(varlist)
      if varlist{1}(1) ~= '-'
         [varlist1, k] = AddFun(varlist{i}, varlist1, k);
      else
         outfun = getoutf(varlist{i});
         if ~isempty(outfun)
            if ~g_grind.statevars.vector
               for h = 1:g_grind.statevars.dim
                  varlist1{k} = i_statevars_names(h);
                  k = k + 1;
               end;
            else
               for i1 = 1:length(g_grind.statevars.vectnames)
                  varlist1{k} = sprintf('outf(''%s'',''%s'')',outfun,g_grind.statevars.vectnames{i1});
                  k = k + 1;
               end;
            end;
         elseif strcmp(varlist{i}, '-extern')
            for h = 1:length(g_grind.externvars)
               varlist1{k} =  g_grind.externvars{h}.name;
               k = k + 1;
            end;
         elseif strcmp(varlist{i}, '-removedata')
            for j = 1:length(g_grind.timevars{No})
               if length(g_grind.timevars{No}{j}) >= 9
                  if strcmpi(g_grind.timevars{No}{j}(1:9), 'observed ')
                     varlist1{k} = g_grind.timevars{No}{j};
                     k = k + 1;
                  end;
               end;
            end;
         elseif strcmp(varlist{i}, '-data')
            if ~isempty(g_data) && ~isempty(g_data.obs)
               for h = 1:g_grind.statevars.dim
                  if ~min(isnan(g_data.obs(:, h)))
                     varlist1{k} =  ['Observed ' i_statevars_names(h)];
                     k = k + 1;
                  end;
               end;
            end;
         elseif strncmpi(varlist{i}, '-def',4)
            if g_grind.statevars.vector
               for h = 1:length(g_grind.statevars.vectnames)
                  varlist1{k} = sprintf('_mean(%s)', g_grind.statevars.vectnames{h});
                  k = k + 1;
               end;
            else
               for h = 1:g_grind.statevars.dim
                  varlist1{k} = i_statevars_names(h);
                  k = k + 1;
               end;
               v1 = getvarlist(No, {'-extern'});
               for h = 1:length(v1)
                  varlist1{k} =  v1{h};
                  k = k + 1;
               end;
               v1 = getvarlist(No, {'-data'});
               for h = 1:length(v1)
                  varlist1{k} =  v1{h};
                  k = k + 1;
               end;
            end;
         else
            [varlist1, k] = AddFun(varlist{i}, varlist1, k);
         end;
      end;
   end;
   if k > 1
      varlist1 = varlist1(1:k - 1);
   else
      varlist1  = {};
   end;
end;

% function [res] = notfound(s)
% global g_grind;
% res = 1;
% for No = 1:length(g_grind.timevars)
%    for i = 1:length(g_grind.timevars{No})
%       if strcmp(g_grind.timevars{No}{i}, s)
%          res = 0;
%          return;
%       end;
%    end;
% end;

function [varlist1, k] = AddFun(fun, varlist1, k)
global g_grind;
if g_grind.statevars.vector
   for h = 1:length(g_grind.statevars.vectnames)
      if strcmp(fun, g_grind.statevars.vectnames{h})
         for j = g_grind.statevars.dims{h}.from:g_grind.statevars.dims{h}.to
            varlist1{k} = i_statevars_names(j);
            k = k + 1;
         end;
         return;
      end;
   end;
end;
varlist1{k} = fun;
k = k + 1;


function addtovars(vars, No)
global g_grind;
if length(g_grind.timevars) < No
   g_grind.timevars{No} = {};
   g_grind.outt{No} = 't';
end;
oudtimevars = g_grind.timevars{No};
k = length(g_grind.timevars{No}) + 1;
for j = 1:length(vars)
   var = vars{j};
   found = 0;
   for i = 1:length(oudtimevars)
      if strcmp(oudtimevars{i}, var)
         found = 1;
         break;
      end;
   end;
   if ~found
      g_grind.timevars{No}{k} =  var;
      k = k + 1;
   end;
end;

function [vars,ntime] = removei(timevars,i, ntime)
if i <= ntime
   if (i > 1) && (i < ntime)
      vars = {timevars{1:i - 1} timevars{i + 1:ntime}};
   elseif i == 1
      vars = timevars(i + 1:ntime);
   elseif i == ntime
      vars = timevars(1:i - 1);
   end;
   ntime = ntime-1;
else
   vars = timevars;
end;

function res = isnumstring(s)
res = 1;
for i = 1:length(s)
   if isempty(strfind('0123456789-', s(i)))
      res = 0;
      return;
   end;
end;

function cleantimevars
global g_grind;
vars = g_grind.timevars;
g_grind.timevars = {};
j = 1;
for i = 1:length(vars)
   if ~isempty(vars{i})
      g_grind.timevars{j} = vars{i};
      j = j + 1;
   end;
end;


function outfun = getoutf(s)
funcs={'_min','_max','_mean','_sum','_median','_std','_cv','_variance','_var','_perc','_shannon'};
outfun = '';
for i = 1:length(funcs)
   fun = funcs{i};
   fun(1) = '-';
   if strncmpi(s, fun, length(fun))
      if strcmpi(fun, '-perc')
         outfun = s(2:end);
      else
         outfun = fun(2:end);
      end;
      return;
   end;
end;


