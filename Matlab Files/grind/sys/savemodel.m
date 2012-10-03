%SAVEMODEL   Save curret parameters and model
%    Save the model and the default parameter settings. After changing a model
%    with the command model this is done automatically.
%
%    Usage:
%    SAVEMODEL- prompts for a filename and saves the model to the 
%    that file.
%    SAVEMODEL FILENAME - saves the model as FILENAME.
%    SAVEMODEL FILENAME 1 - saves the model as FILENAME and overwrites
%    existing files without confirmation.
%  
%    See also savepar, model

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:28 $
function savemodel(afilename, overwrite)
global g_grind;
if nargin < 2
   overwrite = 0;
end;
if ~isfield(g_grind, 'inifile')
   errordlg('No model selected');
   error('GRIND:savemodel:NoModel','No model to save');
end;
if nargin == 0
   if isempty(g_grind.inifile)
      afilename=uiputfile('*.ini', 'Save model as');
   else
      afilename = uiputfile(g_grind.inifile, 'Save model as');
   end;
   overwrite = 1;
end;
overwrite = i_checkstr(overwrite);
%if length(g_grind.model) == 0
%   errordlg('No model specified');
%   error('savemodel');
%end;
if iscell(afilename)
   afilename = afilename{1};
end;
if isempty(strfind(afilename, '.'))
   afilename = [afilename '.ini'];
end;
fid = fopen(afilename, 'r');
if (fid > 0) && ~overwrite
   but=questdlg(['Ok to overwrite ' ,afilename, '?'],'Saving model','Yes','No','Cancel','Yes');
   fclose(fid);
else
   but = 'Yes';
end;
LF = sprintf('\n');
if strcmp('Yes', but)
   [path, name, ext] = fileparts(afilename);
   if exist(afilename, 'file')
      ext = ['.~', ext(2:length(ext))];
      copyfile(afilename, fullfile(path,[name ext]));
   end;
   if ~isempty(strfind(afilename, '.all'))
      saveall(afilename);
      return;
   end;
   fid = fopen(afilename, 'w');
   fwrite(fid, ['%model', LF]);
   for i = 1:length(g_grind.model)
      if ~isempty(g_grind.model{i})
         fwrite(fid, [char(g_grind.model{i}) LF]);
      end;
   end;
   fwrite(fid, ['%commands' LF]);
   for i = 1:length(g_grind.commands)
      if ~isempty(g_grind.commands{i})
         fwrite(fid, [char(g_grind.commands{i}) LF]);
      end;
   end;
   if isfield(g_grind, 'scheme') && ~isempty(g_grind.scheme)
      fwrite(fid, ['%scheme' LF]);
      for i = 1:length(g_grind.scheme)
         if ~isempty(g_grind.scheme{i})
            fwrite(fid, [char(g_grind.scheme{i}) LF]);
         end;
      end;
   end;
   fclose(fid);
   disp(['Model saved to ', afilename]);
   g_grind.inifile = afilename;
elseif strcmp('No', but)
   afilename=uiputfile('*.ini', 'Save model as');
   savemodel(afilename, 1);
end;
function saveall(afilename)
%g_version = 3;
%who global
%eval(i_globalstr(who('global')));
save(afilename, g_grind, '-mat');
disp(['Model and current state saved to ', afilename]);
