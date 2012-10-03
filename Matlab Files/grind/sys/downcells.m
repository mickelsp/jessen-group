%DOWNCELLS   shift down the elements of a matrix/vector
%  Efficient function to shift the elements of a matrix one row down. This
%  function is used to access neighbors of a matrix or vector state variable. 
%  Optionally, the first row is neighboring to the last row. 
%
%  Usage:
%  DOWNCELLS(N,BORDER) - shift the matrix N 1 position down. If BORDER is 1 the first row 
%  borders to itself, if BORDER is 0 the first row borders the last. 
%  DOWNCELLS(N,2,VALUE) - if BORDER is 2 the boundary condition VALUE is used at the border.
% 
%  See also model, rightcells, leftcells, upcells

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:26 $
function Result=downcells(N,bordered,aValue)
if nargin<2
   bordered=1;
end;
if (bordered==2)&&(nargin<3)
   error('GRIND:downcells:ArgError','Downcells: No boundary condition given');
end;
s=size(N,1);
if bordered==1
   Result=[N(1,:);N(1:s-1,:)];
elseif bordered==0
   Result=[N(s,:);N(1:s-1,:)];
elseif bordered==2
   Result=[ones(1,size(N,2))+aValue;N(1:s-1,:)];
end;


