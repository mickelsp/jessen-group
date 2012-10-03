%ONSTART   Function that should be run once before each simulation
%   If some initiation should be done before each simulation, use onstart. 
%   This function is run before each new simulation. Similar things can be
%   done with setevent.
%   
%   Usage:
%   ONSTART FUN - run function FUN before each run (FUN may also be a single line of 
%   commands). (if there was already another function the new function is added)
%   ONSTART -LIST - list the commands.
%   ONSTART -CLEAR - clear all commands.
%   
%   See also setevent, modelpanel

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:28 $
function onstart(varargin)
global g_grind;
if nargin == 0
   for j = 1:length(g_grind.onstart.funs)
      evalin('base', g_grind.onstart.funs{j});
   end;
   return;
end;
fun=strtrim(sprintf('%s ',varargin{:}));
if strncmpi(fun, '-c', 2)
   g_grind.onstart.funs = {};
   return;
elseif strncmpi(fun, '-l', 2)
   if isempty(g_grind.onstart.funs)
      disp('no onstart commands');
   else
      fprintf('onstart %s\n',g_grind.onstart.funs{:})
   end;
   return;
end;
i = length(g_grind.onstart.funs) + 1;
for j = 1:length(g_grind.onstart.funs)
   if strcmp(g_grind.onstart.funs{j}, fun)
      i = j;
      warnging('MATLAB:onstart:alreadyin','Function was already in onstart')
   end;
end;
if isempty(strfind(fun,';'))
    fun=sprintf('%s;',fun);
end;
g_grind.onstart.funs{i} = fun;
evalin('base', g_grind.onstart.funs{i});
