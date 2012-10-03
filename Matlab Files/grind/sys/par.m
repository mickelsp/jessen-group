%PAR   Show parameters and model
%   Show the model equations and the values of parameters
%
%   See also stat

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:28 $
function [p,p2] = par(full, varargin)
global g_grind;
if nargin == 0
   full = '';
elseif ~isempty(full) && strcmp(full, '?')
   e = par('');
   [str,ok]=listdlg('liststring',e,'Name','Parameters',...
      'PromptString','Select parameter','SelectionMode',...
      'single','uh',40,'ffs',10,'fus',10);
   if ok
      f = strfind(e{str}, '=');
      if ~isempty(f)
         parname = strtrim(e{str}(1:f(1) - 1));
         val = strtrim(e{str}(f(1) + 1:length(e{str}) - 1));
         val=inputdlg(['Change parameter ' parname],'Change parameter',1,{val});
         assignin('base', parname, str2num(val{1})); %#ok
      end;
   end;
elseif strncmpi(full, '-e',2) || strncmpi(full, 'edit',5)
   i_parcheck;
   i_pardlg;
   return;
elseif strncmpi(full, '-g',2) || strncmpi(full, 'group',5)
    %Create/list parameter groups
   if nargin == 1
      if isfield(g_grind, 'pargroups')
         pg = sort(unique(g_grind.pargroups));
         maxlen = 0;
         pp = {};
         for i = 1:length(pg)
            if length(pg{i}) > maxlen
               maxlen = length(pg{i});
            end;
         end;
         if length(pg) > 1
            pp = cell(size(pg));
            for i = 1:length(pg)
               pars = g_grind.pars(strcmp(pg{i}, g_grind.pargroups));
               pp{i}=sprintf(['par group %' num2str(maxlen+2) 's  %s'],['''' pg{i} ''''],strtrim(sprintf('%s ',pars{:})));
            end;
         end;
         if nargout >= 1
            p = pg;
            if nargout>1
                p2=pp;
            end;
         else
            fprintf('%s\n', pp{:});
         end
      else
        if nargout >= 1
            p={};
            p2={};
        end;
      end;
      return;
   end;
   defgroup = '-';
   if ~isfield(g_grind, 'pargroups')
      g_grind.pargroups = cell(size(g_grind.pars));
      for i = 1:length(g_grind.pargroups)
         g_grind.pargroups{i} = defgroup;
      end;
   end;
   if nargin==3 && iscell(varargin{2})
      pars = varargin{2};
   elseif nargin == 2
      pars = g_grind.pars(strcmp(defgroup, g_grind.pargroups));
      
   else
      pars = varargin(2:end);
   end;  
   selndx = ones(size(pars));
   pargroup = varargin{1};
   for i = 1:length(pars)
      ndx = strcmp(g_grind.pars, pars{i});
      if ~any(ndx)
         selndx(i) = 0;
         warning('GRIND:par:unknown','Parameter %s is unknown',pars{i});
      else
         g_grind.pargroups{ndx} = pargroup;
      end;
   end;
   return;
elseif strcmp(full, '-maxparlen')
   p = 1;
   if isfield(g_grind, 'model')
      for i = 1:length(g_grind.pars)
         if length(g_grind.pars{i}) > p
            p = length(g_grind.pars{i});
         end;
      end;
   end;
   return;
elseif strncmpi(full, '-d', 2)
   if nargin == 1
      varargin = g_grind.pars;
   end
   pc = cell(1, length(varargin));
   for i = 1:length(varargin)
      p = findparcomm(varargin{i});
      if ~isempty(p)
         f = strfind(p,',');
         if ~isempty(f)
            p = p(1:f(end) - 1);
         end;
      end;
      pc{i} = p;
   end;
   if length(pc) > 1
      p = pc;
   end;
   return;
elseif strncmpi(full, '-setu', 5)
   descr=par('-d',varargin{1});
   comm=sprintf('%%%s,%s',descr,varargin{2});
   replaceparcomm(varargin{1},comm);
   return
elseif strncmpi(full, '-set', 4)
   if isempty(strfind(varargin{2},','))
      comm=sprintf('%%%s,%s',varargin{2},par('-u',varargin{1}));
   else
      comm=['%' varargin{2}];
   end;
   replaceparcomm(varargin{1},comm);
   return;
elseif strncmpi(full, '-u', 2)
   if nargin == 1
      varargin = g_grind.pars;
   end
   pc = cell(1, length(varargin));
   for i = 1:length(varargin)
      p = findparcomm(varargin{i});
      if ~isempty(p)
         f = strfind(p,',');
         if isempty(f)
            p = '';
         else
            p = p(f(end) + 1:end);
         end;
      end;
      pc{i} = p;
   end;
   if length(pc) > 1
      p = pc;
   end;
   return;
end;
if nargout == 0
   % disp(' ');
   if strcmpi(full, 'modelonly')
      disp('Current model');
   elseif strcmpi(full, '-size')
      disp('Sizes of parameters/variables');
      if g_grind.statevars.vector
         dispsiz(g_grind.statevars.vectnames);
      else
         dispsiz(g_grind.statevars.names);
      end;
      dispsiz(g_grind.pars);
      return
   else
      disp('Current model and parameters');
   end;
   if ~isfield(g_grind, 'model') || isempty(g_grind.model)
      disp('No model selected');
   else
      for i = 1:size(g_grind.model, 2)
         disp(char(g_grind.model{i}));
      end;
      if ~isempty(g_grind.scheme)
          disp('<a href="matlab: vismod">scheme</a>');
      end;
   end;
end;
maxparlen = par('-maxparlen');

if isfield(g_grind, 'pars')&&~strcmpi(full, 'modelonly')
   pp = cell(1, size(g_grind.pars, 2));
   for i = 1:size(g_grind.pars, 2)
      ppar = evalin('base', g_grind.pars{i});
      if (size(ppar, 1) > 1) || (size(ppar, 2) > 1)
         s=[g_grind.pars{i} ' = [' ];
         siz = size(ppar);
         sh = 0;
         if ~strcmpi(full, 'full')
            if (siz(1) > 10)
               siz(1) = 10;
               sh = 1;
            end;
            if siz(2) > 10
               siz(2) = 10;
               sh = 1;
            end;
         end;
         for j = 1:siz(1)
            for k = 1:siz(2);
               s = sprintf('%s%g, ',s,ppar(j,k));
            end;
            if j < size(ppar, 1)
               s(length(s) - 1) = ';';
               s = sprintf('%s...\n    ', s);
            end;
         end
         s(length(s) - 1) = ']';
         s(length(s)) = ';';
         if sh
            s =  sprintf('%s\n%%(first 10 elements shown)', s);
         end;
         pp{i} = s;
      elseif ~isempty(ppar)
         pp{i}=sprintf(['%-' num2str(maxparlen) 's = %0.5g;'], g_grind.pars{i}, ppar);
      else
         pp{i}=sprintf(['%-' num2str(maxparlen) 's = [];'], g_grind.pars{i});
      end;
      descr = par('-d', g_grind.pars{i});
      if ~isempty(descr)
         pp{i}=sprintf(['%-' num2str(maxparlen+25-8) 's %%%s,%s'],pp{i},descr,par('-u',g_grind.pars{i}));
      end;
   end;
   if g_grind.solver.hasevents
      disp('Discrete events');
      setevent  '-list';
   end;
   if nargout == 1
      p = pp;
   else
      if isfield(g_grind, 'pars')
         pg=par('groups');
         if isempty(pg)
            fprintf('%s\n',pp{:});
         else
             for j=1:length(pg)
                 fprintf('%%%s:\n',pg{j});
                 ndx=strcmp(pg{j},g_grind.pargroups);
                 fprintf('%s\n',pp{ndx});
             end;
         end;                
      end;
      if ~isempty(g_grind.externvars)
         defextern;
      end;
      i_parcheck(1);
   end;
end;
function dispsiz(v)
for i = 1:size(v, 2)
   siz = evalin('base', sprintf( 'size(%s)',v{i}));
   fprintf('%s = [%dx%d]\n', v{i}, siz);
end;
function [parcomm] = findparcomm(p)
global g_grind;
parcomm = '';
lenp = length(p);
for i = 1:length(g_grind.commands);
   s = strtrim(g_grind.commands{i});
   f1 = strfind(s, '%');
   if ~isempty(f1) && (f1(1) < length(s))
      if strncmpi(p, strtrim(s),lenp)
         if (length(s) > lenp) && ~isempty(strfind(s(lenp + 1:end), '='))
            parcomm = s(f1 + 1:end);
            return;
         end;
      end;
   end;
end;
function replaceparcomm(p,comm)
global g_grind;
lenp = length(p);
for i = 1:length(g_grind.commands);
   s = strtrim(g_grind.commands{i});
   if strncmpi(p, s,lenp)
      f1 = strfind(s, '%');
      if isempty(f1)
          g_grind.commands{i}=sprintf('%s    %s',g_grind.commands{i},comm);
      else
          g_grind.commands{i}=[s(1:f1(1)-1) comm];
      end;
      return;
   end;
end;
