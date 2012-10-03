%SUMNEIGHBORS   Sum all 4 or 8 neighbors
%
%  Usage:
%  SUMNEIGBORS(N,4,BORDER) - Sum 4 neighbors. If BORDER is 1 the first row 
%  borders to itself, if BORDER is 0 the first row borders the last. 
%  SUMNEIGBORS(N,8) - Sum 8 neigbors, BORDER =0 by default. 
%
%  See also model, rightcells, leftcells, upcells, downcells

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:29 $
function Result=sumneighbors(N,nneighbors,bordered)
if nargin<2
   nneighbors=4;
end;
if nargin<3
   bordered=0;
end;
s=size(N);
if bordered
   N1=zeros(s(1),1);
   N2=zeros(1,s(2));
   L=[N1,N(:,1:s(2)-1)];
   R=[N(:,2:s(2)),N1];
   Result=[N(2:s(1),:);N2]+[N2;N(1:s(1)-1,:)]+R+L;
   if nneighbors==8
      Result=Result+[L(2:s(1),:);N2]+[N2;L(1:s(1)-1,:)]+[R(2:s(1),:);N2]+[N2;R(1:s(1)-1,:)];
   end;
else
   L=[N(:,s(2)),N(:,1:s(2)-1)];
   R=[N(:,2:s(2)),N(:,1)];
   Result=[N(2:s(1),:);N(1,:)]+[N(s(1),:);N(1:s(1)-1,:)]+R+L;
   if nneighbors==8
      Result=Result+[L(2:s(1),:);L(1,:)]+[L(s(1),:);L(1:s(1)-1,:)]+[R(2:s(1),:);R(1,:)]+[R(s(1),:);R(1:s(1)-1,:)];
   end;
end;
