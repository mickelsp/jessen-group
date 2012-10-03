%LOADDATA   Load data for external variables or state variables
%   Load a data file with data of external variables or state variables (for 
%   data fitting, see optimpars). The first row of the file should contain 
%   the names of the variables (t=time). If there is no time entered, it is 
%   assumed that the data are equally spaced with time step 1. The format of the
%   ASCI file should be TAB delimited, comma delimited or space delimited.
%
%   Usage:
%   LOADDATA - The user can select a file.
%   LOADDATA FILENAME- Loads the data file FILENAME.
%
%   See also setdata, defextern, optimpars

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:27 $
function loaddata(filename)
global g_grind;
if nargin==0
   if isfield(g_grind,'loaddata')
      filter=[g_grind.loaddata.path filesep '*.csv' pathsep '*.dat' pathsep '*.txt' pathsep '*.ini'];
   else
      filter=['*.csv' pathsep '*.dat' pathsep '*.txt' pathsep '*.ini'];
   end;
   [filename,path]=uigetfile(filter, 'Get data file (delimited or CSV)');
   if filename ~= 0
      filename = [path filename];
   else
      return;
   end;
end;
if strncmpi(filename,'-ini',4)
    filename=g_grind.inifile;
end;
if ~isempty(strfind(filename, ';'))
   filename =filename(1:strfind(filename, ';')-1);
end;
    
[varlist, amatrix] = i_loaddata(filename);
setdata(amatrix,varlist);
