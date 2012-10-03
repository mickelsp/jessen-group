%RIGHTCELLS   Shift right the elements of  a matrix
%  Efficient function to shift the elements of a matrix one column to the right. This
%  function is used to access neighbors of a matrix state variable. 
%  Optionally, the first column is neighboring the last column. 
%
%  Usage:
%  RIGHTCELLS(N,BORDER) - shift the matrix N 1 position to the right. If BORDER is 1 the first column 
%  borders to itself, if BORDER is 0 the first column borders the last. 
%  RIGHTCELLS(N,2,VALUE) - if BORDER is 2 the boundary condition VALUE is used at the border.
%
%  See also model, leftcells, upcells, downcells

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:28 $
function Result=rightcells(N,bordered,aValue)
if nargin<2
   bordered=1;
end;
if (bordered==2)&&(nargin<3)
   error('GRIND:rightcells:NoBoundary','Rightcells: No boundary condition given');
end;
s=size(N,2);
if bordered==1
   Result=[N(:,2:s),N(:,s)];
elseif bordered==0
   Result=[N(:,2:s),N(:,1)];
elseif bordered==2
   Result=[N(:,2:s),ones(size(N,1),1)+aValue];
end;

