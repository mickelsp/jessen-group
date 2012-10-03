%VISMOD   Create a model as a Forrester diagram
%
%
%   Select different elements of a diagram to create a Forrester diagram.
%   Elements are:
%   State variables, name and initial values of state variables.
%   Auxiliary variables, name and equation of auxiliary (or help variables.
%   Such variable contains an equation.
%   Parameters, name and values of parameters.
%   External variables, name and values of exernal variables.
%   External or used-defined functions, name of external function.
%   Continuous flow of substance, a part of a differential equation.
%   Discrete flow of substance, a part of a difference equation.
%   Connectors, information flow between components.
%   Cloud, a source or sink
%
%   Usage:
%   VISMOD - create new model or edit the current model.
%   VISMOD INIFILE - edit inifile.
%   VISMOD -clear - Clears the current model from memory (and deletes currently temporary files).
%
%
%   Note:
%   VISMOD works only on Microsoft Windows systems.
%   
%   
%   See also model, modelpanel

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:29 $
function vismod(amodel)
global g_grind;
if ~strncmpi(computer,'PCWIN',5)
    error('GRIND:vismod:NoPCWIN','"vismod" works only on PC systems, please use <a href="matlab:model">model</a> instead to create a model');    
end;
if (nargin==1)&&(strncmp(amodel,'-c',2))
    model(amodel); % clear the current model
    return;
end;
g=grindpath(2);
if isempty(g_grind)&&~strcmp(pwd, g)
   cd(g);
   fprintf('Changed directory to %s\n',g);
end;
tmpf = 'visualgrind%d.tmp';
i = 0;
tmpfile = sprintf(tmpf, i);
while exist(fullfile(grindpath, tmpfile), 'file')  && (i < 10)
   i = i + 1;
   tmpfile = sprintf(tmpf, i);
end;
if nargin == 0
   if ~isempty(g_grind);
      s = g_grind.inifile;
      amodel = sprintf('curfile#%s', g_grind.inifile);
      savepar(fullfile(grindpath, tmpfile), 1);
      g_grind.inifile = s;
   else
      amodel = 'new';
   end;
else
   if isempty(strfind(amodel, '.'))
      amodel = [amodel '.ini'];
   end;
   [pathname, iscurdir] = findgrindfile(amodel);
   if isempty(pathname)
      errordlg(['Inifile: ' amodel ' does''t exist']);
      error('GRIND:vismod:NoInifile','Inifile "%s" doesn''t exist',model);
   elseif ~iscurdir
      disp(['Could not find the model "' amodel '" in "' cd '"']);
      disp(['Changing to:"' pathname '"']);
      cd(pathname);
   end;
end;
oldpath = cd;
cd(grindpath);
if (length(amodel) > 8) && strcmp(amodel(1:8), 'curfile#')
   name = amodel;
   ext = '';
   p = oldpath;
else
   [p, name, ext] = fileparts(amodel);
end;
if isempty(p)
   p = oldpath;
end;
s = sprintf('visualgrind vismod "%s" "%s" "%s"', p, [name, ext], tmpfile);
%disp(s);
disp('Close the vismod window to continue working in MATLAB');
dos(s);
cd(grindpath);
fid = fopen(tmpfile, 'r');
if fid > 0
   line = myfgetl(fid);
   if strcmp(line, 'OK')
      oldpath = myfgetl(fid);
      inifile = myfgetl(fid);
      fclose(fid);
      delete(tmpfile);
      cd(oldpath);
      use(inifile);
      newpath=fileparts(inifile);
      if ~isempty(newpath)
         cd(newpath);
      end;
   else
      disp('Cancelled')
      fclose(fid);
      delete(tmpfile);
      cd(oldpath);
   end;
else
   error('GRIND:vismod:commerror','Failed to communicate with visualgrind.exe');
end;
function s = myfgetl(fid)
s = fgetl(fid);
if int16(s(end)) == 13
   s = s(1:end - 1);
end;

