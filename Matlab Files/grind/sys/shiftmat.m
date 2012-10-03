%SHIFTMAT - shift the elements of a matrix
%   Shift the elements of a matrix to right, left and/or up and down. This
%   function is used to access neighbors of a matrix or vector state variable. 
%   The first row is neighboring to the last row and the first column is neighboring 
%   the last column. 
%
%  Usage:
%  SHIFTMAT(NX,NY) - shift the matrix NX positions to the right and NY down (if NX is negative, 
%  the matrix shifts to left and if NY is negative it shifts up) 
% 
%  See also: 
%  MODEL
function A1=shiftmat(A,nx,ny)
s=size(A);
if nargin<3
   if s(1)==1
      ny=nx;
      nx=0;
   else
      ny=0;
   end;
end;
if nx<0
   nx=nx+s(2);
end;
if ny<0
   ny=ny+s(1);
end;
if nx>0
   A1=[A(:,s(2)-nx+1:s(2)),A(:,1:s(2)-nx)];
   if ny>0
     A1=[A1(s(1)-ny+1:s(1),:);A1(1:s(1)-ny,:)];
  end;
elseif ny>0
   A1=[A(s(1)-ny+1:s(1),:);A(1:s(1)-ny,:)];
else
   A1=A;
end;
