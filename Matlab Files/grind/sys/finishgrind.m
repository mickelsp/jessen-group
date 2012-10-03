%FINISHGRIND
%   This function should be called automatically if MATLAB is closed.
%   (Add it to a user-defined finish.m function).
%   It cleans-up the CURR_ODE* files.
%
%   Usage:
%   FINISHGRIND = delete the current odefile and set back the counter.
%   FINISHGRIND -reset  = reset the session counter.
function finishgrind(doreset)
global g_grind;
if nargin == 0
   doreset = [];
end;
if isfield(g_grind, 'odefile')
   odefilename=g_grind.odefile;
   gmex('-d');
   if strcmpi(doreset, '-reset')
      if isempty(g_grind)
         setdefaults(1);
      end;
      j = length(odefilename);
      while (j > 0) && ~isempty(intersect(odefilename(j), '1234567890'))
         j = j - 1;
      end;
      odefilename = [odefilename(1:j) '%d.m'];
      for j = 0:100
         f = fullfile(grindpath ,sprintf(odefilename, j));
         if exist(f, 'file')
            delete(f);
         end;
      end;
      odefile = 'curr_pot%d.m';
      for j = 0:100
         f = fullfile(grindpath ,sprintf(odefile, j));
         if exist(f, 'file')
            delete(f);
         end;
      end;
      odefile = 'tmpplt%d.m';
      for j = 0:100
         f = fullfile(grindpath ,sprintf(odefile, j));
         if exist(f, 'file')
            delete(f);
         end;
      end;
   elseif ~isempty(g_grind)
      clear(odefilename)
      f = fullfile(grindpath , odefilename);
      if isempty(strfind(f, '.'))
         f = [f '.m'];
      end;
      if exist(f, 'file')&&~isempty(strfind(odefilename,'curr_ode'))
         delete(f);
      end;
   end;
end;
%rmpath([grindpath filesep 'sys2']);
