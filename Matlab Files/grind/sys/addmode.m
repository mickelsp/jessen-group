%ADDMODE   Sets the simulation mode of the model
%   If the add mode is ON, a new run is appended to the last run, unless the initial
%   conditions change. If add mode is OFF, nothing happens if for instance 
%   TIME is pressed the second time. By default the add mode is OFF.
%
%   Usage:
%   ADDMODE - Toggles add mode ON or OFF. 
%   ADDMODE ON - Set add mode ON. 
%   ADDMODE OFF - Set add mode OFF. 
%   ADDMODE ON -SILENT - doesn't display anything. 
%   ADDMODE -RESET - Resets the current run. 
%
%
%   See also ru, backw, time
%

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:25 $
function addmode(state, opt)
global g_grind g_Y;
i_parcheck;
silent=0;
if nargin==0
   g_grind.solver.addmode=~g_grind.solver.addmode;
else
   if nargin>=2
      if strncmpi(opt,'-s',2)
         silent=1;
      end;
      if strncmpi(opt,'-r',2)
         g_Y=[];
      end;
   end;
   if ~ischar(state)
      state=num2str(state);
   end;
   if strcmpi(state,'on')||strcmpi(state,'1')
      g_grind.solver.addmode=1;
   elseif strcmpi(state,'off')||strcmpi(state,'0')
      g_grind.solver.addmode=0;
   elseif strncmpi(state,'-s',2)
      g_grind.solver.addmode=~g_grind.solver.addmode;
      silent=1;
   elseif strncmpi(state,'-r',2)
         g_Y=[];
         if ~silent
            disp('GRIND Reset');
         end;
         return;
   elseif ~strcmp(state,'?')
      error('GRIND:addmode:UnknownState','?? state must be ON or OFF');
   end;
end;
if ~silent
if g_grind.solver.addmode
   disp('GRIND add mode is on');
else
   disp('GRIND add mode is off');
end;
end;
