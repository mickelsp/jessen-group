function A=drawwaittime(meanperiod,n)
if nargin==1
   n=1;
end;
A=ln(1-rand(n,1))*-meanperiod;
