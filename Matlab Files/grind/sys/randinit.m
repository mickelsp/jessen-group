%RANDINIT   Random initial values
%   Set random (uniformly distributed) initial values.
%   
%   Usage:
%   RANDINIT - chooses random initial points between 0 and 100.
%   RANDINIT MAX - chooses random initial points between 0 and MAX.
%   RANDINIT [MIN MAX] - chooses random initial points between MIN and MAX.
%   RANDINIT VAR - chooses random initial points of state variable VAR between 0 and 100.
%   RANDINIT VAR [MIN MAX] - chooses random initial points of VAR between MIN and MAX.
%   RANDINIT VAR MAX - chooses random initial points of VAR between 0 and MAX.
%
%
%   See also ke, val

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:28 $
function randinit(VarName,range)
global g_grind g_Y;
i_parcheck;
if nargin==0
   range=[0 100];
   VarName=[];
elseif nargin==1
   if ischar(VarName) && isempty(str2num(VarName)) %#ok
      range=[0 100];
   else
      range=i_checkstr(VarName);
      VarName=[];
   end;
elseif nargin==2
   range=i_checkstr(range);
end;
if length(range)==1
   range=[0 range];
end;
if isempty(VarName)
   N0=drawuniform(g_grind.statevars.dim,1,range(1),range(2));
else
   N0=i_initvar;
   [sfrom, sto] = i_statevarnos(VarName);
   if ~isempty(sfrom)
       N0(sfrom:sto)=drawuniform(sto-sfrom+1,1,range(1),range(2));
   end;
end;
i_keep(N0);
g_Y=[];
