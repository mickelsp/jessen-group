%SOLVER   Set the solver and settings.
%   For difference equations there is only one method available.
%   Enter the following information:
%   - Solver (rk4 Euler ode45 ode23 ode115 ode15S ode23S ode23T ode23TB)
%   - settings for the solver, see: MATLAB description or rk4/ euler
%
%   See also rk4, euler, ode45, ode23, ode15i, ode15s, ode23s, ode23t, ode23tb

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:29 $
function res = solver(solv, p1, p2)
global g_grind;
i_parcheck;
     if g_grind.solver.haslag
        solverlist = {'dde23'};
     elseif g_grind.solver.isdiffer
         solverlist = {'i_differ'};
     else
         solverlist = {'rk4','euler','ode45','ode23','ode113','ode15s','ode15i', 'ode23S', 'ode23t', 'ode23tb','ode78','ode87'};
     end;
     wrappers= {'dde23','ddesol';'ode15i','i_ode15sol'};
 if nargout == 1
   if strcmpi(solv, 'step')
      res = g_grind.solver.opt.MaxStep;
      if isempty(res)
         error('GRIND:solver:NoStepSize','No step size defined for the current solver ("%s")\n', g_grind.solver.name);
      end;
      return;
   elseif strcmpi(solv, 'name')
      i=find(strcmpi(wrappers(:,2),g_grind.solver.name));
      if isempty(i)
          res = g_grind.solver.name;
      else
          res=wrappers{i,1};
      end;
      return
   elseif strcmpi(solv, 'list')
      res = solverlist;
      return
   end;
end;
if nargin > 0
   if solv == '?'
      if isempty(g_grind.solver.opt.MaxStep)
         fprintf('solver %s %g %g\n', g_grind.solver.name, g_grind.solver.opt.RelTol, g_grind.solver.opt.AbsTol);
      else
         fprintf('solver %s %g\n', g_grind.solver.name, g_grind.solver.opt.MaxStep);
      end;
      return;
   end;
   n=find(strcmpi(solverlist,solv));
   if ~isempty(n)
 %     if ~(g_grind.solver.isdiffer || g_grind.solver.haslag)
         g_grind.solver.name = solverlist{n};
         i=find(strcmpi(wrappers(:,1),g_grind.solver.name));
         if ~isempty(i)
            g_grind.solver.name =wrappers{i,2};
         end;   
  %   end;
      if (n <= 2) && (nargin > 1)
         g_grind.solver.opt.MaxStep = i_checkstr(p1);
      elseif nargin == 3
         g_grind.solver.opt.RelTol = i_checkstr(p1);
         g_grind.solver.opt.AbsTol = i_checkstr(p2);
      end;
   end;
else
   i_solverdlg;
end

