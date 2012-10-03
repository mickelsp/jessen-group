%UPDATEGRIND   Download newest version from website and install this
% version, overwriting the current version. 
%
%  Usage: 
%  updategrind
%
%  See also setupgrind

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:29 $
function updategrind
curdir=pwd;
try
   cd(grindpath)
   cd ..
   setupgrind('-update');
   cd(curdir)
catch err
   cd(curdir)
   rethrow(err)
end;