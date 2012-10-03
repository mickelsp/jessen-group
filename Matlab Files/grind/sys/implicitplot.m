%IMPLICITPLOT   Implicit plot of an equation
%   Create a two-dimensional implicit plot of an equation, which may
%   include parameters and state variables (no 'functions'). The
%   function is drawn in the same plot that funplot uses.
%   
%   Usage:
%   IMPLICITPLOT EQUATION VARX VARY - create an implicit plot of 
%   VARX versus VARY of EQUATION, which is a function of both VARX and
%   VARY.
%   IMPLICITPLOT EQUATION VARX VARY XRANGE YRANGE - a range for the
%   X-axis and/or the Y-axis may be supplied.
%   
%   Examples:
%   IMPLICITPLOT X^2+Y^2=1 X Y [-1 1] - this function draws a circle.
%   IMPLICITPLOT A*r*(1-A/K)-g*Z*(A/(A+h)) A Z [0.001 10] [0.001 10] draws one of 
%   the nullclines of a algae-zooplankton model.
%   
%   See also funplot, funplot3

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:27 $
function implicitplot(afun, varx, vary, axrange, ayrange, npoints, issurf)
global g_grind;
g_npoints=100;
if nargin<7
   issurf=0; %issurf is used internally by funplot3
end;
if nargin < 3
   if nargin == 0
      afun = '';
   end;
   prompt={'Implicit function','variable on x axis','variable on y axis','range x','range y','number of raster points'};
   if isfield(g_grind, 'implicitplot')
      answer = g_grind.implicitplot;
      if isempty(answer)
         error('GRIND:implicitplot:cancelled','Cancelled');
      elseif ~isempty(afun)
         answer{1} = afun;
      end;
   else
      answer = {afun,'x','y','[0 100]','[0 100]','100'};
   end;
   answer = inputdlg(prompt, 'Implicit function plot', 1, answer);
   if isempty(answer) || isempty(answer{1})
      error('GRIND:implicitplot:NoEquation','No equation entered');
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
   g_grind.implicitplot = answer;
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
fun=afun;
i=strfind(fun, '=');
if ~isempty(i)
   fun=[fun(1:i(1)-1) '- (' fun(i(1)+1:length(fun)) ')'];
end;
if isempty(ayrange)
   ayrange=axrange;
end;
eval(i_globalstr({vary}));
g_loc.comm=sprintf('%s=g_loc.y;',vary);
g_l_ran = ayrange(1):(ayrange(2) - ayrange(1))  / g_npoints:ayrange(2);
eval(sprintf('%s=%g;',vary,g_l_ran(1)));
g_V=zeros(g_npoints+1);
g_X=g_V;
g_Y=g_V;
for i=1:g_npoints+1;
   %  g_loc.y=g_l_ran(i);
   eval(sprintf('%s=%g;',vary,g_l_ran(i)));

%   eval(sprintf('%s=%g;',vary,g_l_ran{i}));
%   g_loc.comm=sprintf('global %s;\n%s=%g;',vary,g_l_ran(i));
   [xx,yy]=i_funcalc(fun,varx,axrange,g_npoints);%,g_loc);
   g_V(:,i)=yy';
   g_X(:,i)=xx';
   g_Y(:,i)=g_l_ran(i);
end;
[ID, new] = i_makefig('funplot');
if new
   hold('on');
end;
plotedit('off');
set(ID,'Name','Function plot');
oldhold = ishold;
hold('on')
if issurf
   surf(g_X,g_Y,g_V);
   shading('flat');
   set(gca,'View',[322.5 ,30]);
 else
   contour(g_X,g_Y,g_V,[0 0],'-b');
end;
if oldhold
   hold('on');
else
   hold('off');
end;
i_plotdefaults;
xlabel(varx);
ylabel(vary);
title(afun);
set(gca, 'YLim', ayrange);
set(gca, 'XLim', axrange);
