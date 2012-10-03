%GRINDPATH   Full path of the grind/sys directory
%   Display the path of the GRIND program files

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:27 $
function p = grindpath(sys)
if nargin==0
   sys=1;
end
if sys==1
   p=fileparts(which('grind.m'));
elseif sys ==0
%  the root of the path of the grind ini files, you may change these lines
  root=fileparts(which('grind.m'));
  root=root(1:length(root)-4);
  grindroot=root(1:length(root)-6);
  p = {root,grindroot};
elseif sys==2
  root=fileparts(which('grind.m'));
  p=root(1:length(root)-4);
end;
