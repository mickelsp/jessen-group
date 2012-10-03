%FUNPLOT3   Plot a 3D surface function
%   Plot some equation in a surface plot (including parameters, state variables, time and
%   functions). For state variables the current initial value is used. 
%
%   Usage:
%   FUNPLOT3 FUN VARX VARY XRANGE YRANGE - FUN is a function which is 
%   plotted on the y-axis ; VARX is the dependent variable for 
%   the X axis; VARY is the dependent variable for 
%   the Y axis; XRANGE is the range of which the variable VAR1 is 
%   varied. YRANGE sets the range of the Y axis. 
%   FUNPLOT3 - without arguments, a dialog box appears.
% 
%   Examples:
%   funplot3 'x^2+y^2' x y [0 10] [0 10]
%   
% 
%   See also funcs, implicitplot, funplot

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:27 $
function funplot3(afun, varx, vary, axrange, ayrange, npoints)
global g_grind;
g_npoints=100;
%if nargin<7
%   issurf=0; %issurf is used internally for funplot3
%end;
if nargin < 3
   if nargin == 0
      afun = '';
   end;
   prompt={'3D function','variable on x axis','variable on y axis','range x','range y','number of raster points'};
   if isfield(g_grind, 'funplot3')
      answer = g_grind.funplot3;
      if isempty(answer)
         error('GRIND:funplot3:Cancelled','Canceled');
      elseif ~isempty(afun)
         answer{1} = afun;
      end;
   else
      answer = {afun,'x','y','[0 100]','[0 100]','100'};
   end;
   answer = inputdlg(prompt, '3D function plot', 1, answer);
   if isempty(answer) || isempty(answer{1})
      error('GRIND:funplot3:NoEquations','no equation entered');
   elseif isempty(answer{2})
      varx = symvar(answer{1});
      varx = varx{1};
   else
      varx = answer{2};
   end;
   if isempty(answer{3})
      vary = symvar(answer{1});
      vary = vary{2};
   else
      vary = answer{3};
   end;
   g_grind.funplot3 = answer;
   afun = answer{1};
   axrange = i_checkstr(answer{4});
   if isempty(axrange)
      axrange = [0.1 100];
   end;
   ayrange = i_checkstr(answer{5});
   g_npoints= i_checkstr(answer{6});
else
   if nargin < 4
      axrange = [0.1 100];
   else
      axrange = i_checkstr(axrange);
   end;
   if nargin < 5
      ayrange = [];
   else
      ayrange = i_checkstr(ayrange);
   end;
   if nargin > 5
      g_npoints = i_checkstr(npoints);
   end;
end;
i=strfind(afun, '=');
if ~isempty(i)
   afun=afun(i(1)+1:length(afun));
end;
implicitplot(afun, varx, vary, axrange, ayrange, g_npoints, 1);

