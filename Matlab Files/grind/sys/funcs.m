%FUNCS   Edit functions
%   Edit/add "functions" (or auxiliary variables) for usage on the axes or for function plots
%   (funplot). The model is not changed. The functions may not have arguments.
%
%   Examples:
%   cons=g*Z*A/(A+H);
%   You may use IF FOR WHILE and SWITCH statements:
%   'Tent function:'
%   if (x(t)<=0.5)
%       fx=r*x(t)
%   else
%       fx=r-r*x(t)
%   end
%   x(t+1)=fx
%   Note that for modelling equations like this, the usage of intermediary variables
%   is required.  
%
%   Usage:
%   FUNCS - opens a screen to edit functions
%   FUNCS FUN1 FUN2 - adds FUN1 and FUN2 to the functions without opening the edit screen.
%
%   See also funplot, ax, out

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:27 $
function funcs(varargin)
global g_grind;
i_parcheck;
eoln = sprintf('\n');
if nargin>0
   if ~iscell(varargin)
      varargin={varargin};
   end;
   for i=1:length(varargin)
      if isempty(strfind(varargin{i},';'))
         g_grind.funcs=[g_grind.funcs eoln varargin{i} ';'];
      else
         g_grind.funcs=[g_grind.funcs eoln varargin{i}];
      end;
   end;
   updatefuncs;
   return;
end;
if isempty(g_grind.funcs)
   g_grind.funcs=sprintf('\n');
end;
answer=inputdlg({'Edit functions'},'Functions for output (doesn''t change model)',20,{g_grind.funcs});
if ~isempty(answer)
   a = answer{1};
   g_grind.funcs = [];
   for i = 1:size(a, 1)
      a1 = strtrim(a(i, :));
      if ~isempty(a1)
         if isempty(strfind(a1,';'))
            g_grind.funcs = [g_grind.funcs a1 ';' eoln];
         else
            g_grind.funcs = [g_grind.funcs a1 eoln];
         end;
      end;
   end;
   g_grind.funcs = g_grind.funcs(1:length(g_grind.funcs) - length(eoln));
   updatefuncs;
end;
disp(g_grind.funcs);
function updatefuncs
global g_grind;
funcs=str2cell(g_grind.funcs);
if isfield(g_grind.funcnames,'dims')
   g_grind.funcnames=rmfield(g_grind.funcnames,'dims');
end;
for i=1:length(funcs)
    f=strfind( funcs{i},'=');
   if ~isempty(f)
      funname = strtrim(funcs{i}(1:f(1) - 1));
      j=1;
      while (j<=length(g_grind.funcnames.names)) &&~strcmp(g_grind.funcnames.names{j},funname)
         j=j+1;
      end;
      if j>length(g_grind.funcnames.names)
         g_grind.funcnames.names=[g_grind.funcnames.names {funname}];
      end;
   end;
end;
clear g_funcs;
