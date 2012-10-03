%COPYFIG   Make a copy of a figure in MATLAB. 
%   Note that when you copy a figure to a figure with a high number, 
%   it is not erased with the command era. Subsequently, the figure can
%   be combined with the command combfig.
%
%   Usage:
%   COPYFIG TONR - Copy the current figure to a new figure with number TONR.
%   COPYFIG FROMNR TONR - Copy figure FROMNR to TONR.
%
%   See also combfig, era

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:26 $
function copyfig(from,to)
if ~exist('i_checkstr','file')
    addpath([grindpath filesep 'sys2']);
end;    
if nargin==1
   to=i_checkstr(from);
   from=gcf;
elseif nargin==2
   from=i_checkstr(from);
   to=i_checkstr(to);
else
   error('GRIND:copyfig:ArgError','No figures specified');
end;
if ishandle(to)
   delete(to);
end;
h=figure(to);
if ishandle(from)
   ax=findobj(from,'type','axes');
   copyobj(ax, h);
else
   error('GRIND:copyfig:UnknownFig','Copyfig: source figure doesn''t exist');
end;

