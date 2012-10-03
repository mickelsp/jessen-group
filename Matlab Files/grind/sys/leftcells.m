%LEFTCELLS   shift left the elements of  a matrix
%  Efficient function to shift the elements of a matrix one column to the left. This
%  function is used to access neighbors of a matrix state variable. 
%  Optionally, the first column is neighboring the last column. 
%
%  Usage:
%  LEFTCELLS(N,BORDER) - shift the matrix N 1 position to the left. If BORDER is 1 the first column 
%  borders to itself, if BORDER is 0 the first column borders the last. 
%  LEFTCELLS(N,2,VALUE) - if BORDER is 2 the boundary condition VALUE is used at the border.
%br>
% 
%  See also model, rightcells, upcells, downcells

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:27 $
function N=leftcells(N,bordered,aValue)
if nargin<2
   bordered=1;
end;
if (bordered==2)&&(nargin<3)
   error('GRIND:leftcells:NoBoundary','Leftcells: No boundary condition given');
end;
s=size(N,2);
% bordered is not periodic
if bordered==1
   N=[N(:,1),N(:,1:s-1)];
elseif bordered==0
   N=[N(:,s),N(:,1:s-1)];
elseif bordered==2 %bordered==2  fixed value
   N=[zeros(size(N,1),1)+aValue,N(:,1:s-1)];
end;

