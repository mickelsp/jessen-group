%FUNPLOT   Plot any 2D equation (including parameters)
%   Plot some equation (including parameters, state variables, time and
%   auxiliary variables). For state variables the current initial value is used. 
%
%   Usage:
%   FUNPLOT FUN VARX XRANGE YRANGE LOGXY- FUN is an equation which is 
%   plotted on the y-axis ; VARX is the dependent variable for 
%   the X axis; XRANGE is the range of which the variable VAR1 is 
%   varied. YRANGE (optional) sets the range of the Y axis. LOGXY 
%   (optional) makes the scale x-axis or the y-axis logarithmic: 
%   values: 'logx' horizontal axis log scale, 'logy' vertical axis 
%   log scale, 'logxy' both log scales.
%   FUNPLOT FUN VAR MIRROR VARRANGE FUNRANGE LOGXY - if MIRROR is 
%   Y then the function is plotted on the x axis and the dependent
%   variable on the y-axis.
%   FUNPLOT VARY FUN YRANGE XRANGE - If the variable and function
%   are exchanged, the function is plotted on the X axis and VARY is
%   the dependent variable for the Y axis. The first range is the 
%   range for the -axis now. 
%   FUNPLOT FUN VARX - if AXRange is omitted, a range of [0.1
%   100] is assumed.
% 
%   Examples:
%   funplot 'A*r*(1-A/K)' A [0 10] ('normal' plot)
%   funplot 'A*r*(1-A/K)' A Y [0 10] ('mirrored' plot)
% 
%   See also funcs, implicitplot, funplot3

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:27 $
function l_HH = funplot(varargin)
%function funplot(afun,g_l_avar,[exchange],[axrange, ayrange],logxy, npoints);
%
global g_grind;
if ~exist('i_use', 'file')
   addpath([grindpath filesep 'sys2']);
end;
%use weird names to avoid that a parameter is changed
g_npoints = 1000;
g_nargin = length(varargin);
defx = 'x';
if g_nargin > 2
   if ischar(varargin{3})
      if ~ismember(upper(varargin{3}), 'NY?')
         varargin = [varargin(1:2), {'N', varargin{3:g_nargin}}];
         g_nargin = g_nargin + 1;
      else
     %    varargin{3} = 'N';
      end;
   else
      varargin = [varargin(1:2), {'N', varargin{3:g_nargin}}];
      g_nargin = g_nargin + 1;
      varargin{3} = 'N';
   end;
   if (g_nargin > 3)
      logxy = i_checklog(varargin{4});
      if ~isempty(logxy)
         varargin = [varargin(1:3), {'[0.01 100]', varargin{4:g_nargin}}];
         g_nargin = g_nargin + 1;
      end;
   else
      varargin=[varargin,{'[0.01 100]','[]',''}];
   end;
   if (g_nargin > 4)
      logxy = i_checklog(varargin{5});
      if ~isempty(logxy)
         varargin = [varargin(1:4), {'[]', varargin{5:g_nargin}}];
         g_nargin = g_nargin + 1;
      end;
   else
      varargin=[varargin,{'[]',''}];
      g_nargin = g_nargin + 2;
   end;
   if (g_nargin < 6)
      varargin = [varargin, {''}];
%      g_nargin = g_nargin + 1;
   end;
else
   if (g_nargin == 1)
      [varargin{1}] = checkeq(varargin{1});
      vars = symvar(varargin{1});
      if (length(vars) == 1);
         varargin = [varargin, vars(1)];
         g_nargin = g_nargin + 1;
      elseif length(vars) > 1
         if ~any(strcmp(vars, defx))
            defx = vars{1};
         end;
      end;
      clear vars;
   end;
   if (g_nargin  < 2)
      if g_nargin == 1
         g_l_afun = varargin{1};
      else
         g_l_afun = '';
      end;
      prompt = {'Function', ...
         'Independent variable','Exchange x and y axis? (y/n)','Range for independent variable','Range for function',...
         'Logarithmic scales (logx/logy/logxy)?'};
      if isfield(g_grind, 'funplot')
         answer = g_grind.funplot;
         if ~isempty(g_l_afun)
            answer{1} = g_l_afun;
         end;
      else
         answer={g_l_afun,defx,'n','[0.1 100]','[]',''};
      end;
      answer = inputdlg(prompt, 'Function plot', 1, answer);
      if isempty(answer) || isempty(answer{1})
         error('GRIND:funplot:NoEquation','No equation entered');
      else
         varargin = answer;
 %        g_nargin = 6;
      end;
      g_grind.funplot = answer;
      clear answer prompt defx;
   else
      varargin={varargin{1:2},'N','[0.01 100]','[]',''};
   end;
end;
g_l_afun = checkeq(varargin{1});
g_l_avar =  varargin{2};
exchange = varargin{3};
axrange = i_checkstr(varargin{4});
ayrange = i_checkstr(varargin{5});
logxy = i_checklog(varargin{6});
clear varargin g_nargin;

%doexchange = 0;
if isempty(exchange)
   [g_l_avar, g_l_afun, doexchange] = checkfun(g_l_avar, g_l_afun, isfunction);
elseif (exchange(1) == '?')
   [g_l_avar, g_l_afun, doexchange] = checkfun(g_l_avar, g_l_afun, isfunction);
else
   doexchange = strcmpi(exchange(1), 'Y');
end;
[g_l_ran, g_l_res] = i_funcalc(g_l_afun, g_l_avar, axrange, g_npoints);
if doexchange
   h1=i_funplot(g_l_res,g_l_ran, g_l_afun,g_l_avar,ayrange,axrange,logxy);
else
   h1=i_funplot(g_l_ran,g_l_res, g_l_avar,g_l_afun,axrange,ayrange,logxy);
end;   
if nargout == 1
   l_HH = h1;
end;
%++++++++++FUNCTIONS++++++++++
function [g_l_avar, g_l_afun, isfunction] = checkfun(g_l_avar, g_l_afun, isfunction)
global g_grind;

%check if there is an operand in g_l_avar
for i = 1:length(g_l_avar)
   if ~isempty(strfind('=*-/+^([.\|&', g_l_avar(i))) %No convert
      isfunction = 1;
      [g_l_avar, g_l_afun] = exchange(g_l_avar, g_l_afun);
      break;
   end
end;
%check if g_l_avar is numeric
isnum = 1;
for i = 1:length(g_l_avar)
   if isempty(strfind('1234567890-.', g_l_avar(i)))
      isnum = 0;
   end;
end;
if isnum
   isfunction  = 1;
   [g_l_avar, g_l_afun] = exchange(g_l_avar, g_l_afun);
end;
%check if g_l_avar is a defined function
if isfield(g_grind, 'funcs')
   funcs = str2cell(g_grind.funcs);
   for i = 1:length(funcs)
      f = strfind(char(funcs{i}), '=');
      if ~isempty(f) && strcmp(strtrim(funcs{i}(1:f(1) - 1)), g_l_avar)
         isfunction = 1;
         [g_l_avar, g_l_afun] = exchange(g_l_avar, g_l_afun);
         break;
      end;
   end;
end;
g_l_afun = checkeq(g_l_afun);
return;

function [g_l_afun] = checkeq(g_l_afun)
f=strfind(g_l_afun, '=');
if ~isempty(f)
   g_l_afun = g_l_afun(f(1) + 1:length(g_l_afun));
end;
return;

function h = i_funplot(x, y, xaxis, yaxis, axrange, ayrange, logxy)
global g_grind;
if isempty(logxy)
   logxy = [0 0];
end;
[ID, new] = i_makefig('funplot');
if new
   hold('on');
end;
plotedit('off');
set(ID,'Name','Function plot');
oldhold = ishold;
hold('on');
if isfield(g_grind, 'pen') && ~isempty(g_grind.pen)
   co = g_grind.pen.color2;
   pe = g_grind.pen.pen;
else
   co = [0 0 1];
   pe = '-';
end;
if length(y) == 1
   y1 = y(1);
   y = x;
   for i = 1:length(y)
      y(i) = y1;
   end;
end
for i = 1:length(x)
   if ~isreal(x(i))
      x(i) = NaN;
   end;
   if ~isreal(y(i))
      y(i) = NaN;
   end;
end;
h = plot(real(x), real(y), pe, 'Color', co);
xlabel(i_disptext(xaxis));
ylabel(i_disptext(yaxis));
if ~isempty(axrange)
   set(gca, 'XLim', axrange);
end;
if isfield(g_grind, 'pen') && ~isempty(g_grind.pen)
   nextpen;
end;
if ~isempty(ayrange)
   set(gca, 'YLim', i_checkstr(ayrange));
end
plotedit('on');
if logxy(1)
   set(gca,'XScale','Log');
else
   set(gca,'XScale','Linear');
end;
if logxy(2)
   set(gca,'YScale','Log');
else
   set(gca,'YScale','Linear');
end;
if ~oldhold
   hold('off');
end;

