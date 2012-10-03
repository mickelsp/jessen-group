%AX   Define axes of phase plane (variable and range) 
%   Define an axis for the phase plane.
%
%   Usage:
%   AX - show current axis settings.
%   AX AXNAME - clear the axis AXNAME (AXNAME can be X, Y or Z). 
%   AX AXNAME [] [LOW HIGH] - set the limits of axis AXNAME.
%   AX AXNAME VAR1 [LOW HIGH] - set the function/variable of
%   AXNAME and set the limits. VAR1 may be a state variable, a 
%   parameter, an auxiliary variable (see funcs), or an MATLAB 
%   expression with these variables.
%   AX AXNAME1<>AXNAME2 - exchange two axes.
%
%   Examples: 
%   ax x A  - sets the state variable 'A' on the x axis.
%   ax y K [0 10] - sets the parameter 'K' on the y axis with a
%   range of 0-10
%   ax z K*A - sets the expression 'K*A' on the z axis.
%   ax z [] [0 100] - sets the range of the z axis to 0-100
%
%   ax z - clears the z axis.
%   ax x<>y - exchanges x and y axis.
% 
%   See also null, null3, phas, ru, paranal, outf

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:26 $
function ax(axname, axvar, axlim, axmax)
global g_grind;
if nargin == 4
   axlim = [i_checkstr(axlim), i_checkstr(axmax)];
end;
if nargin < 3
   axlim = [];
else
   axlim = i_checkstr(axlim);
end;
if nargin < 2
   axvar = ' ';
end;
if (nargin >= 1) && strncmpi(axname, '-def', 4)
   irhs = g_grind.statevars.dim;
   if isempty(g_grind.xaxis.var) && (irhs > 0)
      g_grind.xaxis.var = char(i_statevars_names(1));
   end;
   if isempty(g_grind.yaxis.var) && (irhs > 1)
      g_grind.yaxis.var = char(i_statevars_names(2));
   end;
   if isempty(g_grind.zaxis.var) && (irhs > 2)
      g_grind.zaxis.var = char(i_statevars_names(3));
   end;
   if nargin == 1
      ax('?');
   end;
   return;
end;
if strcmp(axvar, '[]')
   axvar = [];
end;

if (nargin==2)
   f1=strfind(axvar,'[');
   f2=strfind(axvar,']');
   if ~isempty(f1)&&~isempty(f2)
      axlim = i_checkstr(axvar);
      axvar = [];
   end;
end;
if (nargin == 0) || strcmp(axname, '?')
   disp(' ');
   disp ('Current axis settings:');
   if ~isfield(g_grind, 'xaxis')||isempty(g_grind.xaxis) || isempty(g_grind.yaxis) || isempty(g_grind.zaxis)
      disp('currently no axis defined');
   else
      if nargin == 0
         prompt = {'Variable/function on axis x', ...
            'Range axis x', ...
            'Variable/function on axis y', ...
            'Range axis y', ...
            'Variable/function on axis z', ...
            'Range axis z'};
         answer = {char(g_grind.xaxis.var), sprintf('[%g %g]', ...
            g_grind.xaxis.lim) , ...
            char(g_grind.yaxis.var), sprintf('[%g %g]', g_grind.yaxis.lim) , ...
            char(g_grind.zaxis.var), sprintf('[%g %g]', g_grind.zaxis.lim)};
         if isempty(g_grind.yaxis.lim), answer{4} = ''; end;
         if isempty(g_grind.zaxis.lim), answer{6} = ''; end;
         answer = inputdlg(prompt, 'Axes', 1, answer);
         if ~isempty(answer)
            g_grind.xaxis.var = answer{1};
            g_grind.xaxis.lim = i_checkstr(answer{2});
            g_grind.yaxis.var = answer{3};
            g_grind.yaxis.lim = i_checkstr(answer{4});
            g_grind.zaxis.var = answer{5};
            g_grind.zaxis.lim = i_checkstr(answer{6});
         end;
      end;
      printax('x', g_grind.xaxis.var , g_grind.xaxis.lim);
      printax('y', g_grind.yaxis.var , g_grind.yaxis.lim);
      printax('z', g_grind.zaxis.var , g_grind.zaxis.lim);
   end;
else
   switch axname
    case '1'
      axname = 'x';
    case '2'
      axname = 'y';
    case '3'
      axname = 'z';
   end;
   if strcmpi(axname, 'x')
      if ~isempty(axvar)
         g_grind.xaxis.var = strtrim(axvar);
      end;
      if ~isempty(axlim)
         g_grind.xaxis.lim = axlim;
      end;
   elseif strcmpi(axname, 'y')
      if ~isempty(axvar)
         g_grind.yaxis.var = strtrim(axvar);
      end;
      if ~isempty(axlim)
         g_grind.yaxis.lim = axlim;
      end;
   elseif strcmpi(axname, 'z')
      if ~isempty(axvar)
         g_grind.zaxis.var = strtrim(axvar);
      end;
      if ~isempty(axlim)
         g_grind.zaxis.lim = axlim;
      end;
   elseif strcmpi(axname,'x<>y')||strcmpi(axname,'y<>x')
      hlp = g_grind.yaxis;
      g_grind.yaxis = g_grind.xaxis;
      g_grind.xaxis = hlp;
   elseif strcmpi(axname,'x<>z')||strcmpi(axname,'z<>x')
      hlp = g_grind.zaxis;
      g_grind.zaxis = g_grind.xaxis;
      g_grind.xaxis = hlp;
   elseif strcmpi(axname,'y<>z')||strcmpi(axname,'z<>y')
      hlp = g_grind.zaxis;
      g_grind.zaxis = g_grind.yaxis;
      g_grind.yaxis = hlp;
   else
      if isempty(axlim)
         axlim = str2num(axvar); %#ok %not str2double
      end;
      if (~isempty(axlim)) && (length(axlim) > 1)
         if strcmp(axname, g_grind.xaxis.var)
            g_grind.xaxis.lim = axlim;
         elseif strcmp(axname, g_grind.yaxis.var)
            g_grind.yaxis.lim = axlim;
         elseif strcmp(axname, g_grind.zaxis.var)
            g_grind.zaxis.lim = axlim;
         else
            error('GRIND:ax:VarUnknown','Axis "%s" does not exist',axname);
         end;
      end;
   end
end
if isfield(g_grind,'xaxis')
   if g_grind.xaxis.lim(2) < g_grind.xaxis.lim(1)
      h = g_grind.xaxis.lim(2);
      g_grind.xaxis.lim(2) = g_grind.xaxis.lim(1);
      g_grind.xaxis.lim(1) = h;
   end;
   if ~isempty(g_grind.yaxis.lim) && (g_grind.yaxis.lim(2) < g_grind.yaxis.lim(1))
      h = g_grind.yaxis.lim(2);
      g_grind.yaxis.lim(2) = g_grind.yaxis.lim(1);
      g_grind.yaxis.lim(1) = h;
   end;
   if ~isempty(g_grind.zaxis.lim) && (g_grind.zaxis.lim(2) < g_grind.zaxis.lim(1))
      h = g_grind.zaxis.lim(2);
      g_grind.zaxis.lim(2) = g_grind.zaxis.lim(1);
      g_grind.zaxis.lim(1) = h;
   end;
   g_grind.xaxis.var = outf('changeshortcut', g_grind.xaxis.var);
   g_grind.yaxis.var = outf('changeshortcut', g_grind.yaxis.var);
   g_grind.zaxis.var = outf('changeshortcut', g_grind.zaxis.var);
end;
function printax(axname, var, lim)
varno = i_getno(var);
var = outf('changeshortcut', var);
if ~isempty(varno.no)
   fprintf('ax %s %s [%g %g]\n',char(axname), char(var) ,lim)
else
   fprintf('ax %s "%s" is not a valid axis\n',char(axname), char(var));
end;




