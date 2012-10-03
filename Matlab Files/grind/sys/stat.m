%STAT   Display the complete status of GRIND (settings, parameters and state variables)
%   Display the complete status of GRIND (global variables, model, parameters
%   and state variables)
%
%   See also par, ax, val

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:29 $
function stat(full)
global g_grind t g_Y g_t;
if nargin == 0
   full = '';
end;
disp(' ');
if ~isfield(g_grind, 'statevars')
   disp('No model selected');
else
   disp('Global variables of GRIND:');
   if isempty(g_grind.statevars)
      disp('g_grind.statevars     = {}');
   elseif g_grind.statevars.vector
      fprintf('g_grind.statevars       = %s\n', i_cell2str(g_grind.statevars.vectnames));
   else
      fprintf('g_grind.statevars       = %s\n', i_cell2str(g_grind.statevars.names));
   end;
   for No=1:length(g_grind.timevars)
      fprintf('g_grind.timevars{%d}       = %s\n',No, i_cell2str(g_grind.timevars{No}));
   end;
   fprintf('g_grind.pars            = %s\n', i_cell2str(g_grind.pars));
   fprintf('g_grind.ndays           = %0.5g\n', g_grind.ndays);
   fprintf('g_grind.tstep           = %0.5g\n', g_grind.tstep);
   if isempty(g_grind) || isempty(g_grind.pen)
      disp('g_grind.pen          = []');
   else
      fprintf('g_grind.pen.i           = %d\n', g_grind.pen.i);
   end;
   fprintf('g_grind.inifile         = ''%s''\n', g_grind.inifile);
   fprintf('g_grind.odefile         = ''%s''\n' ,g_grind.odefile);
   if isempty(g_grind.solver.name)
      disp('g_grind.solver    = []');
   else
      fprintf('g_grind.slowdown        = %d\n', g_grind.slowdown);
      fprintf('g_grind.truncate        = %d\n', g_grind.truncate);
      fprintf('g_grind.solver.name     = ''%s''\n', g_grind.solver.name);
      fprintf('g_grind.solver.isdiffer = %d\n', g_grind.solver.isdiffer);
      fprintf('g_grind.solver.iters    = %d\n', g_grind.solver.iters);
      fprintf('g_grind.solver.backwards= %d\n', g_grind.solver.backwards);
      fprintf('g_grind.diffto          = %s\n', g_grind.diffto);
   end;
   if ~isempty(g_grind.Jacobian)
      for i = 1:size(g_grind.Jacobian, 1)
         for j = 1:size(g_grind.Jacobian, 2)
            fprintf('g_grind.Jacobian(%d,%d)   = ''%s''\n',i,j, char(g_grind.Jacobian(i,j)));
         end
      end
   end;
   fprintf('t = %0.5g\n', t);
   fprintf('g_Y and g_t = [%dx%d] and [%dx%d] matrices with sim. results\n' , size(g_Y), size(g_t));
   ax ?;
   disp(' ');
   disp('Defined functions:')
   disp(g_grind.funcs);
   par(full);
   val(full);
end;

