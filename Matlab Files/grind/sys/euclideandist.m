%EUCLIDEANDIST   Euclidean distance between vectors
%  Distance between two vectors or between matrixes.
%
%  Usage:
%  D=EUCLIDEANDIST(A,B) - if A and B are vectors 
%  the result is the distance between those vectors
%  if A and B are matrixes, compare the columns of
%  the matrixes.
%  D=EUCLIDEANDIST(A,B,dim) - if dim is 2 then compare
%  the rows of both matrices, if dim is 1 compare columns.

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:26 $
function dist = euclideandist(A, B, dim)
if nargin < 3
   dim = 1;
end;
if dim == 1
   dim2 = 2;
else
   dim2 = 1;
end;
SizeA = size(A);
SizeB = size(B);
if (SizeA(1) == SizeB(1)) && (SizeA(2) == SizeB(2))
     dist = sqrt(sum((A - B).^2,dim));
elseif (SizeA(dim2) == 1) && (SizeB(dim) == SizeA(dim))
   dist = zeros(SizeB(dim2),1);
   if dim==1 
      dist=dist';
   end;
   for i = 1:SizeB(dim2)
      if dim == 1
         dist(i) = sqrt(sum((A - B(:,i)).^2));
      else
         dist(i) = sqrt(sum((A - B(i, :)).^2));
      end;
   end;
elseif (SizeB(dim2) == 1) && (SizeB(dim) == SizeA(dim))
   dist = zeros(SizeA(dim2),1);
   if dim==1 
      dist=dist';
   end;
   for i = 1:SizeA(dim2)
      if dim == 1
         dist(i) = sqrt(sum((A(:,i) - B).^2));
      else
         dist(i) = sqrt(sum((A(i, :) - B).^2));
      end;
   end;
else
   error('GRIND:euclideandist:dimagree','Matrix dimensions must agree');
end
