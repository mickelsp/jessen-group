function dN1=truncN(N,dN)
dN1=dN;
for i=1:20
if N(i)<0.001
%    dN1(i)=-N(i);
dN1(i)=0;
end;
end;
return;
