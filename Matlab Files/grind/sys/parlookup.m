%PARLOOKUP   find a value in the first column of a table
%  This function finds a key value in the first column of a table. It returns the remainder.
%  
%  
%  Usage:
%  PARLOOKUP(Table,Key,extrapolate), lookup a key in the table (first column) with or without
%  extrapolation outside the data range 
%
%  Example:
%  A=[1 3;5 6;10 5;15 3];
%  parlookup(A,4) returns 5.25 (interpolated)
%  parlookup(A,20) returns 0 (outside range)
%  parlookup(A,20,1) returns 3 (extrapolated)
%
%  See also insim

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:28 $
function par = parlookup(tabl, key1, extrapolate, NaNisZero)
if nargin < 4
   NaNisZero = 1;
end;
[N, C] = size(tabl);
if key1 < tabl(1, 1) || key1 > tabl(N, 1)
   if (nargin < 3) || ~extrapolate
      par = NaN .* zeros(1, size(tabl, 2) - 1);
   elseif key1 < tabl(1, 1)
      par = tabl(1, 2:C);
   else
      par = tabl(N, 2:C);
   end;
else
   i = 1;
   while (i<N) && tabl(i, 1) <= key1
      i = i + 1;
   end;
   a = (tabl(i, 2:C) - tabl(i - 1, 2:C)) ./ (tabl(i, 1) - tabl(i - 1, 1));
   b = tabl(i - 1, 2:C) - a .* tabl(i - 1, 1);
   par = a .* key1 + b;
end;
if NaNisZero
   par(isnan(par)) = 0;
end;

