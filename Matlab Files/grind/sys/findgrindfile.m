%FINDGRINDFILE - searches grind subdirectories for file/directory
%
%   Usage:
%   [pathname, iscurdir] = findgrindfile(afile)
function [pathname, iscurdir] = findgrindfile(afile)
pathname = '';
%paths={};
dddir = dir(afile);
if ~isempty(dddir)
   iscurdir = 1;
   if dddir(1).isdir
      pathname=[cd filesep afile];
   else
      pathname=pwd;
   end;
 else
   iscurdir=0;
   grindp = grindpath(0);
   for i = 1:length(grindp)
      pathname = findinsubdirs(grindp{i}, afile);
      if ~isempty(pathname)
            return;
      end;
   end;
end;
function pathnam = findinsubdirs(adir, afile)
dddir = dir([adir filesep afile]);
if ~isempty(dddir)
   if dddir(1).isdir
      pathnam=[adir filesep afile];
   else
      pathnam = adir;
   end;
else
   dddir = dir(adir);
   for i = 3:length(dddir)
      if dddir(i).isdir
         pathnam = findinsubdirs([adir filesep dddir(i).name], afile);
         if ~isempty(pathnam)
            return;
         end;
      end;
   end;
   pathnam = [];
end;
