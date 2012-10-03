% RU   Run the model
%    Run the model and update the currently selected plot (the other plots 
%    that need to be updated are closed). If there is no plot selected, a 
%    1D, 2D, 3D phase plane is opened for 1D, 2D and >2D models respectively.
%    This command can also be used to force a new run even if no parameters 
%    were changed.
%
%    Usage:
%    RU - run with the number of days that is in the global variable g_grind.ndays.
%    Use the command simtime to set the default number of days.
%    RU N - run the model for N days.
%
%    See also simtime, time, phas, null, null3, addmode

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:28 $
function ru(ndays, N0)
global t g_grind;
i_parcheck;
if nargin == 0
   ndays = g_grind.ndays;
else
   ndays = i_checkstr(ndays);
end;
if nargin < 2
   N0 = i_initvar;
end
i = 1;
while ~ishandle(i) && (i < i_figno('maxno'))
   i = i + 1;
end;
if ishandle(i)
   curplot = gcf;
elseif isempty(g_grind.yaxis.var)
   curplot = i_figno('phase1');
   null;
elseif isempty(g_grind.zaxis.var)
   curplot = i_figno('phase2');
   i_makefig('phase2');
   set(gca,'Xlim',g_grind.xaxis.lim);
   set(gca,'Ylim',g_grind.yaxis.lim);
else
   curplot = i_figno('phase3');
   i_makefig('phase3');
   set(gca,'Xlim',g_grind.xaxis.lim);
   set(gca,'Ylim',g_grind.yaxis.lim);
   set(gca,'Zlim',g_grind.zaxis.lim);
   set(gca, 'View', [322.5, 30]);
end;
if (g_grind.drawnow)
   if (curplot == i_figno('phase2'));
      if ~(isempty(i_varno(g_grind.xaxis.var)) || isempty(i_varno(g_grind.yaxis.var)))
         g_grind.solver.opt.OutputFcn = str2func('i_odephas');
      end
   elseif curplot == i_figno('phase3')
      if ~(isempty(i_varno(g_grind.xaxis.var)) || isempty(i_varno(g_grind.yaxis.var)) ...
            ||    isempty(i_varno(g_grind.zaxis.var)))
         g_grind.solver.opt.OutputFcn = str2func('i_odephas');
      end
   end;
else
   g_grind.solver.opt.OutputFcn = str2func('i_odespeed');    
end;
% if g_grind.version.isoctave
%   i_ru(g_grind.odefile, t, ndays, N0, 1);
% else
  oldpointer=get(curplot,'pointer');
  set(curplot, 'pointer', 'watch');
  i_ru(g_grind.odefile, t, ndays, N0, 1);
  set(curplot, 'pointer', oldpointer);
  g_grind.solver.opt.OutputFcn = [];
%end;
i_phas(curplot, 1);
