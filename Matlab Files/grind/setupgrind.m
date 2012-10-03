%SETUPGRIND   Setup GRIND
%   Setup script for GRIND, adds the grind system directory to
%   the default path. This file needs only to be run once after
%   installing GRIND.
%
%
%
%   See also Installing

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:29 $
function setupgrind(option)
url = 'http://content.alterra.wur.nl/webdocs/internet/aew/downloads/';
doanyway = 0;
uninstall = 0;
addstartup = 0;
warning('off','backtrace');
if nargin == 1
   if strncmpi(option, '-d', 2)
      doanyway = 1;
   end;
   if strncmpi(option, '-up', 3)
      updategrindfromweb(url);
   elseif strncmpi(option, '-u', 2)
      uninstall = 1;
   end;
   if strncmpi(option, '-s', 2)
      addstartup = 1;
   end;
end;
k = which('grind.m');
if ~isempty(k) && ~doanyway
   thegrindpath = fileparts(k);
   if uninstall
      rmpath(thegrindpath);
%       try
%          [release, system] = getrelease;
%       catch %#ok
%          release = 11;
%          system = 'matlab';
%       end;
%      isoctave = strcmpi(system, 'octave');
%       if release == 11
%          rmpath(fullfile(thegrindpath, 'R11'))
%       end;
%       if release <= 12
%          rmpath(fullfile(thegrindpath, 'R12'))
%       end;
%       if isoctave
%          rmpath(fullfile(thegrindpath, 'octave'))
%       end;
%      try
         savepath; %adds permanently
%       catch %#ok
%          path2rc;
%       end;
      fprintf('Grind uninstalled from %s (no files removed)\n', thegrindpath);
   else
     if exist_newversion(url)<1
         warning('GRIND:setupgrind:newerversion','There is a newer version available, use <a href="matlab: updategrind">updategrind</a> to update.');
     end;
     fprintf('Grind was already installed in %s\n', thegrindpath);
     if addstartup
         doaddstartup;
      end;
   end;
   return;
elseif uninstall
   fprintf('Grind was not installed\n');
   return;
end;
try
   try
      mroot = OCTAVE_PATH;
   catch %#ok
      mroot = matlabroot;
   end;
catch %#ok
   mroot = '';
end;
thegrindpath=fullfile(mroot,'work','grind','sys');
k = which(fullfile(thegrindpath, 'grind.m'));
if isempty(k)
   thegrindpath = fullfile(pwd, 'sys');
end;
addpath(thegrindpath);
[release, system] = getrelease;
isoctave = strcmpi(system, 'octave');
if release == 11
    error('GRIND:setupgrind:R11','MATLAB version R11 is not supported, download the R11 version of GRIND');
%   addpath(fullfile(thegrindpath, 'R11'));
end;
if release == 12
    error('GRIND:setupgrind:R12','MATLAB version R12 is not supported, download the R11 version of GRIND');
%   addpath(fullfile(thegrindpath, 'R12'))
end;
if isoctave
    error('GRIND:setupgrind:octave','Octave is not supported, download the R11 version');
%   addpath(fullfile(thegrindpath, 'octave'));
end;
try
   savepath; %adds permanently
catch %#ok
   path2rc;
end;
w = 0;
if ~canwriteto(grindpath)
   warning('GRIND:setupgrind:writepermission','No permission to write to %s %s %s\n\n',grindpath,...
   'If you have permission in another directory it is possible to run GRIND.',...
   'Preferrably set write permission to the full GRIND directory.');
   w = 1;
end;
p=fullfile(mroot,'toolbox','local');
if ~canwriteto(p)
   warning('GRIND:setupgrind:writepermission','The path cannot be stored as there is no write permission in %s. %s\n\n',...
       p, 'Each session GRIND needs to be reinstalled (or the permission should be changed).');
   w = 1;
end;
if addstartup
   doaddstartup;
end;
if exist_newversion(url)<1
    warning('GRIND:setupgrind:newerversion','There is a newer version available, use <a href="matlab: updategrind">updategrind</a> to update.');
end;
if w
   disp('GRIND installed, but there are some problems (see above)');
else
   disp('Grind installed successfully');
end;

function doaddstartup
thegrindpath = grindpath;
u = userpath;
if length(u) < 4  %R11 gives not the workdirectory
   u = fullfile(matlabroot, 'work');
end;
f = strfind(u, '; ');
if ~isempty(f)
   u = u(1:f(1) - 1);
end;
u = fullfile(u, 'startup.m');
if exist(u, 'file') == 2
%   lines = {};
   fid = fopen(u, 'r');
   hassetup = 0;
   while ~feof(fid)
      line = fgetl(fid);
%      lines = {lines, line};
      if ~isempty(strfind(line, 'setupgrind'))
         hassetup = 1;
      end;
   end;
   fclose(fid);
   if ~hassetup
      fid = fopen(u, 'a');
      fprintf(fid,'\ncd ''%s''\nsetupgrind\n',thegrindpath(1:end-4));
      fclose(fid);
      disp('setupgrind added to startup.m');
   else
      disp('setupgrind was already in startup.m');
   end;
else
   fid = fopen(u, 'w');
   fprintf(fid,'\ncd ''%s''\nsetupgrind\n',thegrindpath(1:end-4));
   disp('startup.m created');
   fclose(fid);
end;
function ok = canwriteto(path)
ok = 0;
oldpath = pwd;
try
   cd(path);
   fid=fopen('t.tmp','wt');
   ok =fid >= 0;
   if ok
      fclose(fid);
   end;
   delete('t.tmp');
   cd(oldpath);
catch %#ok
   cd(oldpath);
end;

function updategrindfromweb(url)
if ~exist('urlread','file')
   error('GRIND:updategrind:MATLABversion','This function works only for newer versions of MATLAB');
end;
[nonewer, web, current] = exist_newversion(url);
if isempty(web)
   error('GRIND:updategrind:noconnection','Cannot download GRIND, is internet connected?')
end;
if nonewer
   warning('GRIND:updategrind:nonewversion','There is no newer version of GRIND available');
else
   fprintf('Removing current version (%s)...',current.date)
   cd(grindpath);
   delete('*.m');
   delete('*.exe');
   cd('sys2')
   delete('*.m');
   dir  '*.m'
   cd ..
   cd ..
   cd ..
   fprintf('Downloading new version (%s)...',web.date)
   unzip([url 'grind.zip']);
   disp('Successfully updated');
end;
fprintf('GRIND version %s\nRelease date: %s\n', web.version, web.date);

function [existnew, web, current] = exist_newversion(url)
existnew = 0;
current = [];
web = [];
try
   h = str2cell(urlread([url 'grindversion.txt']));
catch %#ok
   return;
end;
web.version = h{1};
web.date = h{2};
fid=fopen('use.m','r');
while ~feof(fid)
   line = fgetl(fid);
   f = strfind(line, 'Revision:');
   if ~isempty(f)
      f1 = strfind(line, '$');
      current.version = strtrim(line(f + 9:f1(1) - 1));
      f = strfind(line, 'Date:');
      current.date = strtrim(line(f + 5:f1(end) - 1));
   end;
end;
fclose(fid);
existnew = datenum(web.date)-datenum(current.date) < 1E-4;


