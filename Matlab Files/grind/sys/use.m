%USE   Use model
%   Opens an inifile and use the model. If an wildcard (asterisk) is used,
%   a list of possibilities is displayed
%   
%   Usage:
%   USE opens the model dialog box (similar as MODEL).
%   USE AMODEL - Opens AMODEL if in current directory, else it searches for the file.
%   USE AMOD* - finds files with the mask AMOD*.
%
%   See also model

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:29 $
function use(varargin)
%amodel = [];
%if the model name has spaces it can be read as separate arguments
if isempty(varargin)
   model
   return;
elseif  (nargin==1)&&strncmpi(varargin{1},'-c',2)
   model -clear;
   return;
end;

amodel=strtrim(sprintf('%s ',varargin{:}));
i_startup;
if ~exist('i_use','file')
  addpath([grindpath filesep 'sys2']);
end;
i_use(amodel,1,1,1,1);
par('modelonly');
