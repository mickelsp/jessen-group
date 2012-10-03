%function res=i_rednoise(n,labda,T0,beta,res0);
%n= number of steps, labda,T0,beta = red noise parameters
%res0 (optional) = first value (default res0=T0);
%
function res=i_rednoise(n,labda,T0,beta,res0)
res=repmat(T0,n, 1);
if nargin==5
   res(1)=res0;
end;
for i = 2:n
   res(i) = (1 - 1 / labda) * (res(i - 1) - T0) + T0 + beta * randn;
end
