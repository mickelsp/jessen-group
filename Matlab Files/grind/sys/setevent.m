%SETEVENT   Insert a discrete event in the queue
%    This command is to create a discrete event. A continuous model can be 
%    interrupted by a event (=function) which is scheduled in an event 
%    queue. It is not necessary that this command is done in the model definition, 
%    but it can only be placed in the lower parameter panel. The settings are restored.
%    Example of an event:
%    
%    %Example of an event function
%    function nextt=runevent(t);
%    global A;
%    A=A+4;
%    nextt=t+365;
%
%    Set the result of the function (nextt) to NaN if the event should be deleted.
%
%    Usage:
%    SETEVENT EVENT FIRSTT - EVENT is the name of the event function, FIRSTT is the 
%    first time that the event should take place. FISTT can be a parameter, to get
%    it updated before each run enter a string with the name of the parameter.
%    SETEVENT -CLEAR - Clear all events.
%    SETEVENT -LIST - List the current events.
%    SETEVENT('simpleevent',FIRSTT,COMM,NEXTT) - create a simple event (e.g. change of 
%    state variable or parameter), COMM = string with commands (e.g. 'A=A-1; if A<0.1,
%    A=0.1;end;') NEXTT = time lag between the events.
%
%    See also model
%

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:29 $
function setevent(event, dt,varargin)
if nargin == 0
   i_seteventdlg;
   i_setevent filldlg;
   i_setevent ListButtonDwn;
   return;
elseif nargin==1
   if strncmpi(event,'-c',2)
      i_setevent clear;
      return;
   elseif strncmpi(event,'-l',2)
      i_setevent list;
      return;
   elseif strncmpi(event,'-e',2)
      i_setevent('enable',1);
      return;
   elseif strncmpi(event,'-d',2)
      i_setevent('enable',0);
      return;
   end;
   dt=0;
end;
i_setevent('addnew',event,dt,varargin);

