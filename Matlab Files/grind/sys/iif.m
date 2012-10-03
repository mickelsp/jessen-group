%IIF   Immediate if
%   Conditional function. 
%
%   Usage:
%   IFF(COND,IFTRUE,IFFALSE) if cond is true then the result 
%   is IFTRUE else the result is IFFALSE. All arguments may 
%   contain matrices.
% 
%   Example:
%   x=[0:100]
%   A=iif(x<50,x/50,1) is a short notation of:
%
%   for i=0:100
%     if x(i)<50
%        A(i)=x(i)/50;
%     else
%        A(i)=1;
%     end
%   end;

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:27 $
function res = iif(cond, A1, A2)
if nargin < 3
   error('GRIND:iif:ArgError','"iff" needs three arguments: iif(condition,iftrue,iffalse)');
end
res=ones(size(cond));
if (length(cond) > 1)
   if length(A1)==1
       res(cond)=A1;
   else
       res(cond) = A1(cond);
   end
   if length(A2)==1
       res(~cond)=A2;
   else
       res(~cond) = A2(~cond);
   end;
else 
   if cond
       res = A1;
   else
       res = A2;
   end; 
end; 




