%KE   Keep the final state after a run
%   Fix the final state of the last simulation as new starting
%   point.
%   
%   Usage:
%   KE - saves the final state
%   KE -start (-s) - saves the start state
%   KE T - saves the state closest to time T (no interpolation is used)
%
%   
%   See also val, ru

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:27 $

%   See also:
%   VAL, RU
function ke(opt)
global g_Y g_t g_grind;   
if (nargin==1)
   if strncmpi(opt,'-s',2)
      ndx=1;
   else
      opt=i_checkstr(opt);
      if ~ischar(opt)
         aa=abs(g_t-opt);
         ndx=find(aa==min(aa));
      end;
   end;
else
   ndx=size(g_Y,1);
end;
if ~isempty(g_Y)
   i_keep(g_Y(ndx, :));
   if ~isempty(g_grind.permanent)
      defpermanent('-s',defpermanent('-p',g_t(ndx))); 
   end;
else
   warning('GRIND:ke:norun','No simulation run to keep, use <a href="matlab:time">time</a> or <a href="matlab:ru">ru</a> to simulate first');
end;
