function y = fitalinedangit(coeff,time)
m=coeff(1);
b=coeff(2);

t=time(:);

y = m.*t+b;