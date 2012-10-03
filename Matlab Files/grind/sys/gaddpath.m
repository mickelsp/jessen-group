%gaddpath - simple way to add a directory to the search path
%
%gaddpath stats -sub - include subdirectories
%gaddpath - use dialog box to
function gaddpath(nam, opt)
if nargin < 2
   opt = '-sub';
end;
oldpath = cd;
try
   if nargin == 0
      gcd
   else
      gcd(nam);
   end;
   adir = cd;
   addpath(adir);
   fprintf('<a href="matlab: cd(''%s'')">%s</a> added to the path\n',adir,adir);
   if strncmp(opt, '-s', 2)
       d=dir;
       for i=1:length(d)
           if d(i).isdir && ~strncmp(d(i).name,'...', length(d(i).name))
               adir1=fullfile(adir,d(i).name);
               fprintf('<a href="matlab: cd(''%s'')">%s</a> added to the path\n',adir1,adir1);
               addpath(adir1);
           end;
       end;
   end;
   cd(oldpath);
catch %#ok
   cd(oldpath);
end;
