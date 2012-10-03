function r=smod(i1,i2)
if getrelease==11 %extreme slow mod in R11
   d=i1/i2;
   r=round((d-floor(d))*i2);
else
   r=mod(i1,i2);
end
%d=i1/i2;
%r=round((d-floor(d))*i2);
